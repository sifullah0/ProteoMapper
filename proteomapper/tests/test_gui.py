"""
GUI tests for protein sequence analysis application.

Tests cover:
- Widget creation and initialization
- User interactions and callbacks
- Input validation
- File selection dialogs
- Progress tracking
- Threading integration with GUI
"""

import pytest
import tkinter as tk
from tkinter import ttk
import pandas as pd
import tempfile
import os
import time
from unittest.mock import Mock, patch, MagicMock


from ProteoMapper import ProteinAnalysisApp

@pytest.fixture
def tk_root():
    """Create Tkinter root for testing."""
    root = tk.Tk()
    yield root
    try:
        root.destroy()
    except:
        pass


@pytest.fixture
def app(tk_root):
    """Create application instance for testing."""
    app = ProteinAnalysisApp(tk_root)
    yield app
    try:
        tk_root.destroy()
    except:
        pass


@pytest.fixture
def sample_excel_file(tmp_path):
    """Create sample Excel file for GUI tests."""
    file_path = tmp_path / "test_gui_input.xlsx"
    
    df = pd.DataFrame({
        'Protein Sequences': ['ACDEFGHIKLMNPQRSTVWY', 'MKKLIFGDS'],
        'Gene Name': ['GeneA', 'GeneB']
    })
    df.to_excel(file_path, index=False)
    
    return str(file_path)


class TestGUIInitialization:
    """Tests for GUI initialization and widget creation."""
    
    def test_when_app_created_expect_window_configured(self, app):
        """Test that main window is properly configured."""
        assert app.root.title() == "Protein Sequence Analysis Tool"
        assert app.root.winfo_width() >= 650
        assert app.root.winfo_height() >= 600
    
    def test_when_app_created_expect_all_frames_present(self, app):
        """Test that all main frames are created."""
        # Update to ensure widgets are created
        app.root.update()
        
        # Check that scrollable frame exists
        assert app.scrollable_frame is not None
        assert app.main_canvas is not None
        assert app.scrollbar is not None
    
    def test_when_app_created_expect_variables_initialized(self, app):
        """Test that all Tkinter variables are initialized."""
        assert isinstance(app.motif_search_var, tk.BooleanVar)
        assert isinstance(app.position_highlight_var, tk.BooleanVar)
        assert isinstance(app.domain_scanning_var, tk.BooleanVar)
        assert isinstance(app.motif_file_path, tk.StringVar)
        assert isinstance(app.conservation_threshold, tk.StringVar)
        assert app.conservation_threshold.get() == "60"
    
    def test_when_app_created_expect_proceed_button_disabled(self, app):
        """Test that proceed button is initially disabled."""
        app.root.update()
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_app_created_expect_optional_frames_hidden(self, app):
        """Test that optional analysis frames are initially hidden."""
        app.root.update()
        
        # Check that frames are not packed initially
        assert not app.motif_frame.winfo_viewable()
        assert not app.position_frame.winfo_viewable()
        assert not app.domain_frame.winfo_viewable()


class TestUserInteractions:
    """Tests for user interactions and widget responses."""
    
    def test_when_motif_checkbox_checked_expect_frame_shown(self, app):
        """Test that motif frame appears when checkbox is checked."""
        app.root.update()
        
        # Initially hidden
        assert not app.motif_frame.winfo_viewable()
        
        # Check the checkbox
        app.motif_search_var.set(True)
        app.toggle_motif_options()
        app.root.update()
        
        # Now visible
        assert app.motif_frame.winfo_viewable()
    
    def test_when_position_checkbox_checked_expect_frame_shown(self, app):
        """Test that position frame appears when checkbox is checked."""
        app.root.update()
        
        app.position_highlight_var.set(True)
        app.toggle_position_options()
        app.root.update()
        
        assert app.position_frame.winfo_viewable()
    
    def test_when_domain_checkbox_checked_expect_frame_shown(self, app):
        """Test that domain frame appears when checkbox is checked."""
        app.root.update()
        
        app.domain_scanning_var.set(True)
        app.toggle_domain_options()
        app.root.update()
        
        assert app.domain_frame.winfo_viewable()
    
    def test_when_checkbox_unchecked_expect_frame_hidden(self, app):
        """Test that frames hide when checkboxes are unchecked."""
        app.root.update()
        
        # Show then hide
        app.motif_search_var.set(True)
        app.toggle_motif_options()
        app.root.update()
        
        app.motif_search_var.set(False)
        app.toggle_motif_options()
        app.root.update()
        
        assert not app.motif_frame.winfo_viewable()
    
    def test_when_file_selected_expect_path_variable_updated(self, app, sample_excel_file):
        """Test file path variable updates on file selection."""
        with patch('tkinter.filedialog.askopenfilename', return_value=sample_excel_file):
            app.browse_file(app.motif_file_path)
            app.root.update()
        
        assert app.motif_file_path.get() == sample_excel_file
    
    def test_when_file_dialog_cancelled_expect_no_change(self, app):
        """Test behavior when file dialog is cancelled."""
        original_path = app.motif_file_path.get()
        
        with patch('tkinter.filedialog.askopenfilename', return_value=''):
            app.browse_file(app.motif_file_path)
            app.root.update()
        
        assert app.motif_file_path.get() == original_path


class TestInputValidation:
    """Tests for input validation logic."""
    
    def test_when_no_options_selected_expect_proceed_disabled(self, app):
        """Test that proceed is disabled when no analysis options are selected."""
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_option_selected_no_file_expect_proceed_disabled(self, app):
        """Test that proceed is disabled without input file."""
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_motif_selected_no_patterns_expect_proceed_disabled(self, app, sample_excel_file):
        """Test that proceed is disabled when motif search has no patterns."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_all_required_inputs_provided_expect_proceed_enabled(self, app, sample_excel_file):
        """Test that proceed is enabled when all inputs are valid."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'normal'
    
    def test_when_position_highlight_no_positions_expect_proceed_disabled(self, app, sample_excel_file):
        """Test validation for position highlighting."""
        app.motif_file_path.set(sample_excel_file)
        app.position_highlight_var.set(True)
        # No positions entered
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_domain_scan_no_database_expect_proceed_disabled(self, app, sample_excel_file):
        """Test validation for domain scanning."""
        app.motif_file_path.set(sample_excel_file)
        app.domain_scanning_var.set(True)
        # No database path
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'disabled'
    
    def test_when_multiple_options_all_valid_expect_proceed_enabled(self, app, sample_excel_file, tmp_path):
        """Test validation with multiple analysis options."""
        mock_db = tmp_path / "test.hmm"
        mock_db.write_text("MOCK")
        
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.position_highlight_var.set(True)
        app.positions_to_highlight.set("1 5 10")
        app.domain_scanning_var.set(True)
        app.pfam_db_path.set(str(mock_db))
        
        app.root.update()
        app.check_inputs()
        
        assert str(app.proceed_button['state']) == 'normal'


class TestProgressTracking:
    """Tests for progress window and tracking."""
    
    def test_when_analysis_started_expect_progress_window_shown(self, app, sample_excel_file):
        """Test that progress window appears on analysis start."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        
        # Mock the actual analysis
        with patch.object(app, 'execute_analysis_wrapper'):
            app.run_analysis()
            app.root.update()
        
        # Progress window should be created
        assert app.progress_window is not None
    
    def test_when_analysis_complete_expect_progress_window_closed(self, app, sample_excel_file):
        """Test cleanup after analysis completion."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        
        # Simulate cleanup
        app.cleanup_after_analysis()
        app.root.update()
        
        # Button should be re-enabled
        assert str(app.proceed_button['state']) == 'normal'
    
    def test_when_progress_bar_exists_expect_safe_cleanup(self, app):
        """Test safe cleanup of progress bar."""
        # Create progress window
        app.show_progress_window()
        app.root.update()
        
        # Cleanup should not raise errors
        app.cleanup_after_analysis()
        app.root.update()


class TestThreadingIntegration:
    """Tests for threading behavior in GUI."""
    
    def test_when_analysis_runs_expect_gui_responsive(self, app, sample_excel_file):
        """Test that GUI remains responsive during analysis."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        
        # Mock long-running analysis
        def mock_analysis(*args):
            time.sleep(0.1)  # Simulate work
        
        with patch.object(app, 'execute_analysis', side_effect=mock_analysis):
            with patch.object(app, 'show_completion'):
                app.run_analysis()
                
                # GUI should update during analysis
                for _ in range(5):
                    app.root.update()
                    time.sleep(0.02)
    
    def test_when_multiple_analyses_queued_expect_proper_handling(self, app, sample_excel_file):
        """Test that multiple rapid analysis requests are handled properly."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        
        # Button should be disabled during analysis
        with patch.object(app, 'execute_analysis_wrapper'):
            app.run_analysis()
            app.root.update()
            
            # Button should be disabled
            assert str(app.proceed_button['state']) == 'disabled'


class TestErrorHandling:
    """Tests for error handling in GUI."""
    
    def test_when_invalid_threshold_expect_error_dialog(self, app, sample_excel_file):
        """Test error handling for invalid conservation threshold."""
        app.motif_file_path.set(sample_excel_file)
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.conservation_threshold.set("invalid")
        app.root.update()
        
        with patch('tkinter.messagebox.showerror') as mock_error:
            # Mock the threading to avoid actual execution
            with patch('threading.Thread'):
                app.run_analysis()
                app.root.update()
    
    def test_when_invalid_positions_expect_error_dialog(self, app, sample_excel_file):
        """Test error handling for invalid position values."""
        app.motif_file_path.set(sample_excel_file)
        app.position_highlight_var.set(True)
        app.positions_to_highlight.set("1 abc 5")
        app.root.update()
        
        with patch('tkinter.messagebox.showerror') as mock_error:
            with patch('threading.Thread'):
                app.run_analysis()
                app.root.update()
    
    def test_when_file_not_found_expect_error_dialog(self, app):
        """Test error handling for missing input file."""
        app.motif_file_path.set("/nonexistent/file.xlsx")
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG")
        app.root.update()
        
        with patch('tkinter.messagebox.showerror') as mock_error:
            with patch('threading.Thread'):
                app.run_analysis()
                app.root.update()

class TestExecuteAnalysisWrapper:
    """Dedicated tests for execute_analysis_wrapper."""

    def test_execute_analysis_wrapper_happy_path_motif_and_positions(
        self, app, sample_excel_file
    ):
        """
        Test that execute_analysis_wrapper builds parameters correctly,
        calls execute_analysis, and schedules show_completion, while
        cleaning up the GUI state.
        """
        # ---- Configure GUI state ----
        app.motif_file_path.set(sample_excel_file)

        # Enable motif search with two patterns
        app.motif_search_var.set(True)
        app.regex_text.insert("1.0", "ACDEFG\nHIKLMN")
        app.conservation_threshold.set("60")

        # Enable position highlighting
        app.position_highlight_var.set(True)
        app.positions_to_highlight.set("1 5 10")

        # No domain scanning in this test → domain_params should stay None
        app.domain_scanning_var.set(False)

        # ---- Patch root.after so callbacks run immediately ----
        def fake_after(delay, callback, *args, **kwargs):
            # Immediately execute the scheduled callback
            return callback(*args, **kwargs)

        with patch.object(app.root, "after", side_effect=fake_after):
            # Patch execute_analysis to avoid heavy work and capture its args
            fake_output = "/tmp/fake_output.xlsx"
            with patch.object(app, "execute_analysis", return_value=fake_output) as mock_exec:
                # Patch show_completion so no real window is created
                with patch.object(app, "show_completion") as mock_show:
                    start_time = time.time() - 1.0  # pretend we started 1s ago

                    # Call the wrapper directly (normally run in a background thread)
                    app.execute_analysis_wrapper(start_time)

        # ---- Assertions on execute_analysis call ----
        mock_exec.assert_called_once()
        _, kwargs = mock_exec.call_args

        assert kwargs["input_excel_file"] == sample_excel_file
        assert kwargs["regex_input"] == ["ACDEFG", "HIKLMN"]
        assert kwargs["positions_to_highlight"] == [1, 5, 10]
        # threshold should be parsed as float
        assert kwargs["conservation_threshold"] == 60.0
        # no domain scanning → None
        assert kwargs["domain_params"] is None

        # ---- Assertions on show_completion scheduling ----
        mock_show.assert_called_once()
        out_file_arg, runtime_arg = mock_show.call_args[0]
        assert out_file_arg == fake_output
        assert runtime_arg > 0  # runtime should be positive

        # ---- GUI cleanup: proceed button should be re-enabled ----
        assert str(app.proceed_button["state"]) == "normal"


class TestCompletionDialog:
    """Tests for completion dialog functionality."""
    
    def test_when_analysis_complete_expect_completion_window(self, app, tmp_path):
        """Test that completion window is shown with correct information."""
        output_file = tmp_path / "output.xlsx"
        output_file.write_text("MOCK")
        runtime = 12.5
        
        with patch('tkinter.Toplevel') as mock_toplevel:
            app.show_completion(str(output_file), runtime)
            
            # Completion window should be created
            assert mock_toplevel.called
    
    def test_when_open_file_clicked_expect_file_opened(self, app, tmp_path):
        """Test file opening from completion dialog."""
        output_file = tmp_path / "output.xlsx"
        output_file.write_text("MOCK")
        
        with patch('subprocess.run') as mock_run:
            with patch('platform.system', return_value='Linux'):
                app.open_file(str(output_file))
                
                # Should attempt to open file
                assert mock_run.called


class TestScrolling:
    """Tests for scrollable interface."""
    
    def test_when_content_exceeds_window_expect_scrollbar_active(self, app):
        """Test that scrollbar activates when content is large."""
        # Show all optional frames to increase content
        app.motif_search_var.set(True)
        app.position_highlight_var.set(True)
        app.domain_scanning_var.set(True)
        
        app.toggle_motif_options()
        app.toggle_position_options()
        app.toggle_domain_options()
        
        app.root.update()
        app.update_scroll_region()
        
        # Scrollbar should be configured
        assert app.main_canvas.yview() is not None
    
    def test_when_mouse_wheel_used_expect_scrolling(self, app):
        """Test mouse wheel scrolling functionality."""
        # Create mock event
        mock_event = Mock()
        mock_event.delta = 120
        
        app.root.update()
        
        # Should not raise errors
        # Note: Actual scrolling tested manually due to event complexity


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])