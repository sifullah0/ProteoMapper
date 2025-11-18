import pytest
import pandas as pd
import tempfile
import os
import sys
import shutil
import subprocess

from pathlib import Path

# Directory containing tests/
THIS_DIR = Path(__file__).resolve().parent
# Project root (one level up from tests/)
ROOT_DIR = THIS_DIR.parent

# Ensure project root is on sys.path so `import ProteoMapper` works
root_str = str(ROOT_DIR)
if root_str not in sys.path:
    sys.path.insert(0, root_str)



def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "performance: marks performance/benchmark tests"
    )
    config.addinivalue_line(
        "markers", "hmmer: marks tests requiring HMMER installation"
    )
    config.addinivalue_line(
        "markers", "gui: marks GUI-specific tests"
    )
    config.addinivalue_line(
        "markers", "integration: marks integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks unit tests"
    )
    config.addinivalue_line(
        "markers", "threading: marks threading/concurrency tests"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test items during collection."""
    # Skip HMMER tests if HMMER is not installed
    skip_hmmer = pytest.mark.skip(reason="HMMER not installed")
    
    try:
        subprocess.run(["hmmscan", "-h"], stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE, check=True)
        hmmer_available = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        hmmer_available = False
    
    for item in items:
        if "hmmer" in item.keywords and not hmmer_available:
            item.add_marker(skip_hmmer)


@pytest.fixture(scope="session")
def test_data_dir(tmp_path_factory):
    """Create a temporary directory for test data."""
    data_dir = tmp_path_factory.mktemp("test_data")
    yield data_dir
    # Cleanup
    shutil.rmtree(data_dir, ignore_errors=True)


@pytest.fixture
def sample_sequences():
    """Provide sample protein sequences for testing."""
    return [
        'ACDEFGHIKLMNPQRSTVWY',
        'MKKLIFGDSMQRSTVWY',
        'ACDEFGHIKLMNPQRSTVWY',
        'QRSTACDEFGHIKLMN'
    ]


@pytest.fixture
def sample_gene_names():
    """Provide sample gene names for testing."""
    return ['GeneA', 'GeneB', 'GeneC', 'GeneD']


@pytest.fixture
def sample_dataframe(sample_sequences, sample_gene_names):
    """Create a sample DataFrame for testing."""
    return pd.DataFrame({
        'protein sequences': sample_sequences,
        'gene name': sample_gene_names
    })


@pytest.fixture
def simple_excel_file(tmp_path, sample_sequences, sample_gene_names):
    """Create a simple Excel file for testing."""
    file_path = tmp_path / "simple_test.xlsx"
    
    df = pd.DataFrame({
        'Protein Sequences': sample_sequences,
        'Gene Name': sample_gene_names
    })
    
    df.to_excel(file_path, index=False)
    return str(file_path)


@pytest.fixture
def excel_with_metadata(tmp_path):
    """Create Excel file with metadata rows before header."""
    file_path = tmp_path / "with_metadata.xlsx"
    
    with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
        # Write metadata and header
        metadata_df = pd.DataFrame([
            ['Study Information', ''],
            ['Date:', '2024-01-01'],
            ['', ''],
            ['Protein Sequences', 'Gene Name'],
            ['ACDEFGHIKLMNPQRSTVWY', 'GeneA'],
            ['MKKLIFGDSMQRSTVWY', 'GeneB']
        ])
        metadata_df.to_excel(writer, index=False, header=False)
    
    return str(file_path)


@pytest.fixture
def excel_with_gaps(tmp_path):
    """Create Excel file with sequences containing gaps."""
    file_path = tmp_path / "with_gaps.xlsx"
    
    df = pd.DataFrame({
        'Protein Sequences': [
            'AC--DEFG',
            'HIK-LMN-OP',
            'QR---ST'
        ],
        'Gene Name': ['GeneA', 'GeneB', 'GeneC']
    })
    
    df.to_excel(file_path, index=False)
    return str(file_path)


@pytest.fixture
def large_test_dataset(tmp_path):
    """Create a larger dataset for performance testing."""
    file_path = tmp_path / "large_test.xlsx"
    
    sequences = ['ACDEFGHIKLMNPQRSTVWY'] * 100
    gene_names = [f'Gene{i:03d}' for i in range(100)]
    
    df = pd.DataFrame({
        'Protein Sequences': sequences,
        'Gene Name': gene_names
    })
    
    df.to_excel(file_path, index=False)
    return str(file_path)


@pytest.fixture
def pattern_file(tmp_path):
    """Create a regex patterns file for testing."""
    file_path = tmp_path / "patterns.txt"
    
    patterns = [
        "ACDEFG::Motif1",
        "HIK.*MN::Motif2",
        "[AC]{3,}::PolyAC"
    ]
    
    file_path.write_text('\n'.join(patterns))
    return str(file_path)


@pytest.fixture
def positions_file(tmp_path):
    """Create a positions file for testing."""
    file_path = tmp_path / "positions.txt"
    file_path.write_text("1 5 10 15 20")
    return str(file_path)


@pytest.fixture
def mock_hmm_database(tmp_path):
    """Create a minimal mock HMM database."""
    hmm_file = tmp_path / "mock_pfam.hmm"
    
    hmm_content = """HMMER3/f [3.3.2 | Nov 2020]
NAME  Test_Domain
ACC   PF00001
DESC  Test protein domain
LENG  50
//
"""
    hmm_file.write_text(hmm_content)
    return str(hmm_file)


@pytest.fixture(autouse=True)
def cleanup_after_test(request):
    """Cleanup resources after each test."""
    yield
    # Perform cleanup
    # This runs after each test
    import gc
    gc.collect()


@pytest.fixture(scope="function")
def isolated_directory(tmp_path):
    """Provide an isolated directory for tests that create files."""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    
    original_dir = os.getcwd()
    os.chdir(work_dir)
    
    yield work_dir
    
    os.chdir(original_dir)


# Performance benchmark fixtures
@pytest.fixture
def benchmark_config():
    """Configuration for performance benchmarks."""
    return {
        'rounds': 5,
        'warmup_rounds': 1,
        'iterations': 1
    }


# GUI testing fixtures
@pytest.fixture
def tk_root():
    """Create Tkinter root for GUI tests."""
    import tkinter as tk
    
    root = tk.Tk()
    root.withdraw()  # Hide window during tests
    
    yield root
    
    try:
        root.destroy()
    except:
        pass


# Marker-based fixtures
def pytest_runtest_setup(item):
    """Setup for specific test markers."""
    # Add custom setup for marked tests
    if "slow" in item.keywords:
        # Could add timeout or skip logic here
        pass
    
    if "performance" in item.keywords:
        # Could set up performance monitoring
        pass


def pytest_runtest_teardown(item):
    """Teardown for specific test markers."""
    # Add custom teardown
    pass


# Custom assertion helpers
class Helpers:
    """Helper methods for tests."""
    
    @staticmethod
    def assert_dataframe_equal(df1, df2, check_dtype=True):
        """Assert two DataFrames are equal."""
        pd.testing.assert_frame_equal(df1, df2, check_dtype=check_dtype)
    
    @staticmethod
    def assert_file_exists(file_path):
        """Assert file exists."""
        assert os.path.exists(file_path), f"File not found: {file_path}"
    
    @staticmethod
    def assert_excel_readable(file_path):
        """Assert Excel file is readable."""
        try:
            pd.read_excel(file_path)
        except Exception as e:
            pytest.fail(f"Cannot read Excel file {file_path}: {e}")


@pytest.fixture
def helpers():
    """Provide helper methods to tests."""
    return Helpers()

@pytest.fixture(scope="session")
def toy_hmm():
    return Path(__file__).parent / "fixtures" / "toy.hmm"
