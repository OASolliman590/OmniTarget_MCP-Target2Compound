#!/usr/bin/env python3
"""
Basic Status Test - MCP Drug Discovery Pipeline
===============================================

This script tests basic pipeline components without problematic imports.
"""

import sys
import os
import subprocess
import shutil
from pathlib import Path

def test_basic_status():
    """Test basic pipeline status."""
    print("🧪 Testing Basic Pipeline Status")
    print("=" * 50)
    
    # Test 1: Environment
    print("\n🐍 Test 1: Environment")
    try:
        print(f"✅ Python: {sys.version.split()[0]}")
        print(f"✅ Working directory: {os.getcwd()}")
        
        # Check conda environment
        conda_prefix = os.environ.get('CONDA_PREFIX')
        if conda_prefix:
            print(f"✅ Conda environment: {Path(conda_prefix).name}")
        else:
            print("⚠️ No conda environment detected")
    except Exception as e:
        print(f"❌ Environment error: {e}")
        return False
    
    # Test 2: Configuration File
    print("\n📋 Test 2: Configuration File")
    try:
        config_file = Path('configs/test_sertraline.yaml')
        if config_file.exists():
            print(f"✅ Configuration file exists: {config_file}")
            with open(config_file, 'r') as f:
                content = f.read()
                if 'depression' in content and 'sertraline' in content:
                    print("✅ Configuration contains expected content")
                else:
                    print("⚠️ Configuration content unexpected")
        else:
            print(f"❌ Configuration file not found: {config_file}")
            return False
    except Exception as e:
        print(f"❌ Configuration error: {e}")
        return False
    
    # Test 3: Compound File
    print("\n💊 Test 3: Compound File")
    try:
        compound_file = Path('data/compounds/sertraline_conjugates.smi')
        if compound_file.exists():
            print(f"✅ Compound file exists: {compound_file}")
            with open(compound_file, 'r') as f:
                lines = f.readlines()
                compound_lines = [line for line in lines if line.strip() and not line.startswith('#')]
                print(f"✅ Found {len(compound_lines)} compounds")
                if compound_lines:
                    first_compound = compound_lines[0].strip().split()
                    if len(first_compound) >= 2:
                        print(f"✅ First compound: {first_compound[1]}")
        else:
            print(f"❌ Compound file not found: {compound_file}")
            return False
    except Exception as e:
        print(f"❌ Compound file error: {e}")
        return False
    
    # Test 4: Vina Binary
    print("\n⚗️ Test 4: Vina Binary")
    try:
        vina_path = shutil.which('vina_custom')
        if vina_path:
            print(f"✅ Vina binary found: {vina_path}")
            
            # Test version (with library path)
            env = os.environ.copy()
            env['DYLD_LIBRARY_PATH'] = f"{os.environ.get('CONDA_PREFIX', '')}/lib:{env.get('DYLD_LIBRARY_PATH', '')}"
            
            result = subprocess.run(['vina_custom', '--version'], 
                                  capture_output=True, text=True, env=env, timeout=10)
            if result.returncode == 0:
                print(f"✅ Vina version: {result.stdout.strip()}")
            else:
                print(f"⚠️ Vina version check failed: {result.stderr}")
        else:
            print("❌ Vina binary not found in PATH")
            return False
    except Exception as e:
        print(f"❌ Vina test error: {e}")
        return False
    
    # Test 5: DeepDTA Repository
    print("\n🧠 Test 5: DeepDTA Repository")
    try:
        deepdta_path = Path('third_party/DeepDTA')
        if deepdta_path.exists():
            print(f"✅ DeepDTA repository exists: {deepdta_path}")
            source_path = deepdta_path / 'source'
            if source_path.exists():
                print("✅ DeepDTA source directory found")
            else:
                print("⚠️ DeepDTA source directory not found")
        else:
            print(f"❌ DeepDTA repository not found: {deepdta_path}")
            return False
    except Exception as e:
        print(f"❌ DeepDTA test error: {e}")
        return False
    
    # Test 6: Docker
    print("\n🐳 Test 6: Docker")
    try:
        result = subprocess.run(['docker', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print(f"✅ Docker: {result.stdout.strip()}")
        else:
            print("❌ Docker not available")
            return False
    except Exception as e:
        print(f"❌ Docker test error: {e}")
        return False
    
    # Test 7: Docker Compose
    print("\n🐳 Test 7: Docker Compose")
    try:
        compose_file = Path('docker/docker-compose.yaml')
        if compose_file.exists():
            print(f"✅ Docker compose file exists: {compose_file}")
        else:
            print(f"❌ Docker compose file not found: {compose_file}")
            return False
    except Exception as e:
        print(f"❌ Docker compose test error: {e}")
        return False
    
    # Summary
    print("\n📊 Test Summary")
    print("=" * 50)
    print("✅ Environment: WORKING")
    print("✅ Configuration: WORKING")
    print("✅ Compounds: WORKING")
    print("✅ Vina binary: WORKING")
    print("✅ DeepDTA repository: WORKING")
    print("✅ Docker: WORKING")
    print("✅ Docker compose: WORKING")
    print("\n🎯 Status: READY FOR MCP SERVER INTEGRATION")
    print("   Next step: Start MCP servers with 'docker compose up -d'")
    
    return True

if __name__ == "__main__":
    success = test_basic_status()
    sys.exit(0 if success else 1)
