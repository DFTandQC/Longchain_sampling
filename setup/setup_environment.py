#!/usr/bin/env python3
"""
自动环境配置脚本 - Molecular Cluster Sampling System
Automatic Environment Setup Script

在新计算机上运行此脚本以自动配置环境：
Run this script on a new computer to automatically configure the environment:

    python setup_environment.py

功能 (Features):
1. 检查 Python 版本 (Check Python version)
2. 检查并安装依赖包 (Check and install dependencies)
3. 验证项目结构 (Verify project structure)
4. 运行快速测试 (Run quick smoke test)
5. 显示快速开始指南 (Show quick start guide)
"""

import sys
import subprocess
import os
from pathlib import Path

# ANSI colors for terminal output
class Colors:
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'
    
    @staticmethod
    def supports_color():
        """Check if terminal supports color"""
        return hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

def print_colored(text, color=''):
    """Print colored text if terminal supports it"""
    if Colors.supports_color():
        print(f"{color}{text}{Colors.END}")
    else:
        print(text)

def print_header(text):
    """Print section header"""
    print_colored(f"\n{'='*60}", Colors.BOLD)
    print_colored(f"  {text}", Colors.BOLD + Colors.BLUE)
    print_colored(f"{'='*60}", Colors.BOLD)

def print_success(text):
    """Print success message"""
    print_colored(f"✓ {text}", Colors.GREEN)

def print_warning(text):
    """Print warning message"""
    print_colored(f"⚠ {text}", Colors.YELLOW)

def print_error(text):
    """Print error message"""
    print_colored(f"✗ {text}", Colors.RED)

def check_python_version():
    """检查 Python 版本是否满足要求"""
    print_header("Step 1: Checking Python Version")
    
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    print(f"Current Python version: {version_str}")
    
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print_error(f"Python 3.8+ required, but found {version_str}")
        print("Please upgrade Python: https://www.python.org/downloads/")
        return False
    
    print_success(f"Python {version_str} meets requirements (3.8+)")
    return True

def check_pip():
    """检查 pip 是否可用"""
    print_header("Step 2: Checking pip")
    
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        print_success(f"pip is available: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError:
        print_error("pip is not available")
        print("Please install pip: https://pip.pypa.io/en/stable/installation/")
        return False

def install_dependencies():
    """安装依赖包"""
    print_header("Step 3: Installing Dependencies")
    
    requirements_file = Path(__file__).parent / "requirements.txt"
    
    if not requirements_file.exists():
        print_warning("requirements.txt not found, skipping dependency installation")
        return True
    
    print(f"Installing from: {requirements_file}")
    
    try:
        subprocess.run(
            [sys.executable, "-m", "pip", "install", "-r", str(requirements_file)],
            check=True
        )
        print_success("All dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to install dependencies: {e}")
        return False

def verify_dependencies():
    """验证关键依赖包"""
    print_header("Step 4: Verifying Dependencies")
    
    dependencies = {
        'numpy': '1.20.0'
    }
    
    all_ok = True
    for package, min_version in dependencies.items():
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            print_success(f"{package}: {version} (min: {min_version})")
        except ImportError:
            print_error(f"{package}: NOT INSTALLED (required: {min_version}+)")
            all_ok = False
    
    return all_ok

def verify_project_structure():
    """验证项目结构"""
    print_header("Step 5: Verifying Project Structure")
    
    project_root = Path(__file__).parent
    required_items = {
        'files': ['main.py', 'config.json'],
        'dirs': ['lib', 'monomer', 'docs', 'test']
    }
    
    all_ok = True
    
    # Check required files
    for file in required_items['files']:
        path = project_root / file
        if path.exists():
            print_success(f"Found: {file}")
        else:
            print_error(f"Missing: {file}")
            all_ok = False
    
    # Check required directories
    for dir_name in required_items['dirs']:
        path = project_root / dir_name
        if path.is_dir():
            # Count files in directory
            file_count = len(list(path.glob('*')))
            print_success(f"Found: {dir_name}/ ({file_count} items)")
        else:
            print_error(f"Missing: {dir_name}/")
            all_ok = False
    
    # Create output directory if missing
    out_dir = project_root / 'out'
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
        print_warning("Created missing directory: out/")
    
    return all_ok

def run_smoke_test():
    """运行快速测试"""
    print_header("Step 6: Running Smoke Test")
    
    project_root = Path(__file__).parent
    smoke_test = project_root / 'test' / 'smoke_test.py'
    
    if not smoke_test.exists():
        print_warning("Smoke test not found, skipping")
        return True
    
    print("Running quick validation test...")
    
    try:
        result = subprocess.run(
            [sys.executable, str(smoke_test)],
            capture_output=True,
            text=True,
            timeout=30,
            cwd=str(project_root)
        )
        
        if result.returncode == 0 and 'SMOKE_TEST: OK' in result.stdout:
            print_success("Smoke test PASSED")
            print("\nTest output:")
            print(result.stdout)
            return True
        else:
            print_error("Smoke test FAILED")
            if result.stdout:
                print("\nStdout:")
                print(result.stdout)
            if result.stderr:
                print("\nStderr:")
                print(result.stderr)
            return False
    except subprocess.TimeoutExpired:
        print_error("Smoke test timed out (>30s)")
        return False
    except Exception as e:
        print_error(f"Failed to run smoke test: {e}")
        return False

def show_quick_start():
    """显示快速开始指南"""
    print_header("Setup Complete! Quick Start Guide")
    
    print("""
🚀 Environment setup successful! Try these commands:

1. List available molecules:
   python main.py --list-monomers

2. Generate clusters with default settings:
   python main.py --nseeds 100

3. Select specific molecules:
   python main.py --use "PT,H2SO4,NO2" --nseeds 50

4. Custom molecule counts:
   python main.py --use "PT,H2SO4,NO,NO2" \\
     --counts "PT=2,H2SO4=1,NO=2,NO2=1" --nseeds 200

5. Enable filtering:
   python main.py --use "PT,H2SO4" \\
     --enable_filter --min_contacts 20 --max_rg 25.0 --nseeds 100

6. Run tests:
   python test/tests.py
   python test/test_pt_core.py
   python test/test_filtering.py

📚 Documentation:
   - Full guide: docs/FULL_DOCUMENTATION.md
   - Math details: docs/ALGORITHM.md
   - README: README.md

📁 Output location: out/seeds/

Happy sampling! 🧪
""")

def main():
    """主函数"""
    print_colored("\n" + "="*60, Colors.BOLD)
    print_colored("  Molecular Cluster Sampling System", Colors.BOLD + Colors.BLUE)
    print_colored("  Environment Setup", Colors.BOLD + Colors.BLUE)
    print_colored("="*60 + "\n", Colors.BOLD)
    
    # Run setup steps
    steps = [
        ("Python Version", check_python_version),
        ("pip", check_pip),
        ("Dependencies Installation", install_dependencies),
        ("Dependencies Verification", verify_dependencies),
        ("Project Structure", verify_project_structure),
        ("Smoke Test", run_smoke_test)
    ]
    
    failed_steps = []
    
    for step_name, step_func in steps:
        try:
            success = step_func()
            if not success:
                failed_steps.append(step_name)
        except Exception as e:
            print_error(f"Unexpected error in {step_name}: {e}")
            failed_steps.append(step_name)
    
    # Summary
    print_header("Setup Summary")
    
    if not failed_steps:
        print_success("All setup steps completed successfully!")
        show_quick_start()
        return 0
    else:
        print_error(f"Setup completed with {len(failed_steps)} issue(s):")
        for step in failed_steps:
            print(f"  • {step}")
        print("\nPlease resolve the issues above and run setup again.")
        print("For help, see: docs/FULL_DOCUMENTATION.md")
        return 1

if __name__ == "__main__":
    sys.exit(main())
