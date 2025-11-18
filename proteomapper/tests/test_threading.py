"""
Threading and concurrency tests for protein sequence analysis.
"""

import pytest
import threading
import pandas as pd
from queue import Queue
import sys
import os


from ProteoMapper import (
    clean_and_process_sequences,
    prepare_final_dataframe,
    setup_excel_styles,
    apply_pattern_highlights
)


@pytest.fixture
def sample_dataframe():
    """Create sample DataFrame for threading tests."""
    return pd.DataFrame({
        'protein sequences': ['ACDEFGHIKLMNPQRSTVWY'] * 100,
        'gene name': [f'Gene{i}' for i in range(100)]
    })


class TestThreadSafety:
    """Tests for thread safety of core functions."""
    
    def test_concurrent_sequence_cleaning(self, sample_dataframe):
        """Test that sequence cleaning is thread-safe."""
        results = []
        errors = []
        
        def clean_sequences(df_copy):
            try:
                result = clean_and_process_sequences(df_copy)
                results.append(result)
            except Exception as e:
                errors.append(e)
        
        # Create multiple threads
        threads = []
        for i in range(5):
            df_copy = sample_dataframe.copy()
            thread = threading.Thread(target=clean_sequences, args=(df_copy,))
            threads.append(thread)
            thread.start()
        
        # Wait for all threads
        for thread in threads:
            thread.join(timeout=10)
        
        # All threads should complete without errors
        assert len(errors) == 0
        assert len(results) == 5
        
        # All results should be identical
        for result in results:
            assert len(result) == 100
            assert 'protein sequences_cleaned' in result.columns
    
    def test_concurrent_dataframe_preparation(self, sample_dataframe):
        """Test thread safety of DataFrame preparation."""
        # Clean sequences first
        df_cleaned = clean_and_process_sequences(sample_dataframe.copy())
        
        results = []
        errors = []
        
        def prepare_df(df_copy):
            try:
                result, leading_cols = prepare_final_dataframe(df_copy)
                results.append((result, leading_cols))
            except Exception as e:
                errors.append(e)
        
        threads = []
        for i in range(5):
            df_copy = df_cleaned.copy()
            thread = threading.Thread(target=prepare_df, args=(df_copy,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join(timeout=10)
        
        assert len(errors) == 0
        assert len(results) == 5
        
        # Verify consistency
        first_result, first_cols = results[0]
        for result, cols in results[1:]:
            assert len(result) == len(first_result)
            assert cols == first_cols
    
    def test_no_race_conditions_in_pattern_matching(self):
        """Test for race conditions in pattern matching operations."""
        from openpyxl import Workbook
        
        # Create test data
        df = pd.DataFrame({
            'protein sequences': ['ACDEFGHIKLMNPQRSTVWY'] * 10,
            'gene name': [f'Gene{i}' for i in range(10)]
        })
        df = clean_and_process_sequences(df)
        df_final, leading_cols = prepare_final_dataframe(df)
        
        results = []
        errors = []
        
        def match_patterns():
            try:
                wb = Workbook()
                ws = wb.active
                patterns = [('ACDEFG', 'Motif1')]
                styles = setup_excel_styles()
                
                result = apply_pattern_highlights(ws, df, patterns, leading_cols, styles)
                results.append(result)
            except Exception as e:
                errors.append(e)
        
        threads = []
        for i in range(3):
            thread = threading.Thread(target=match_patterns)
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join(timeout=15)
        
        # Should complete without errors
        assert len(errors) == 0


class TestDataRaceConditions:
    """Tests for potential data race conditions."""
    
    def test_concurrent_dataframe_access(self):
        """Test concurrent access to shared DataFrame."""
        df = pd.DataFrame({
            'protein sequences': ['ACDEFGHIKLMNPQRSTVWY'] * 100,
            'gene name': [f'Gene{i}' for i in range(100)]
        })
        
        results = Queue()
        errors = Queue()
        
        def process_dataframe(thread_id):
            try:
                # Each thread gets its own copy
                df_copy = df.copy()
                df_clean = clean_and_process_sequences(df_copy)
                df_final, _ = prepare_final_dataframe(df_clean)
                results.put((thread_id, len(df_final)))
            except Exception as e:
                errors.put((thread_id, e))
        
        threads = []
        for i in range(10):
            thread = threading.Thread(target=process_dataframe, args=(i,))
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join(timeout=10)
        
        # No errors should occur
        assert errors.qsize() == 0
        assert results.qsize() == 10
        
        # All results should be consistent
        while not results.empty():
            thread_id, length = results.get()
            assert length == 100
    
    def test_no_shared_state_corruption(self):
        """Test that concurrent operations don't corrupt shared state."""
        results = []
        
        def get_styles():
            styles = setup_excel_styles()
            results.append(len(styles))
        
        threads = []
        for i in range(20):
            thread = threading.Thread(target=get_styles)
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join(timeout=5)
        
        # All threads should get same number of styles
        assert all(r == results[0] for r in results)


class TestMemoryInThreads:
    """Tests for memory management in threaded environment."""
    
    def test_memory_cleanup_after_threads(self, sample_dataframe):
        """Test that memory is properly released after threads complete."""
        import gc
        
        def process_data():
            df = sample_dataframe.copy()
            df = clean_and_process_sequences(df)
            df_final, _ = prepare_final_dataframe(df)
            # Data should be garbage collected when thread ends
        
        threads = []
        for i in range(10):
            thread = threading.Thread(target=process_data)
            threads.append(thread)
            thread.start()
        
        for thread in threads:
            thread.join(timeout=10)
        
        # Force garbage collection
        gc.collect()
        
        # All threads should have completed and cleaned up
        active_threads = threading.active_count()
        # Should only be main thread (and possibly daemon threads)
        assert active_threads <= 2


if __name__ == "__main__":

    pytest.main([__file__, "-v"])
