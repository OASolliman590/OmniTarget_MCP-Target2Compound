#!/usr/bin/env python3
"""
Test Current Status - MCP Drug Discovery Pipeline
=================================================

This script tests the current status of the pipeline components
without requiring MCP servers to be running.
"""

import sys
import os
import asyncio
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from orchestrator.schemas.config import RunConfig
from orchestrator.schemas.compounds import Compound
from orchestrator.adapters.deepdta_adapter import DeepDTAAdapter
import yaml

async def test_current_status():
    """Test current pipeline status."""
    print("🧪 Testing Current Pipeline Status")
    print("=" * 50)
    
    # Test 1: Configuration Loading
    print("\n📋 Test 1: Configuration Loading")
    try:
        with open('configs/test_sertraline.yaml', 'r') as f:
            config_data = yaml.safe_load(f)
        config = RunConfig(**config_data)
        print(f"✅ Configuration loaded successfully")
        print(f"   Disease terms: {config.disease_terms}")
        print(f"   Organism: {config.organism}")
        print(f"   Max compounds: {config.compounds.max_compounds}")
    except Exception as e:
        print(f"❌ Configuration error: {e}")
        return False
    
    # Test 2: Compound Loading
    print("\n💊 Test 2: Compound Loading")
    try:
        compounds = []
        with open('data/compounds/sertraline_conjugates.smi', 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 1:
                        smiles = parts[0]
                        name = parts[1] if len(parts) > 1 else f'compound_{line_num}'
                        compound = Compound(name=name, smiles=smiles, compound_id=name)
                        compounds.append(compound)
        
        print(f"✅ Loaded {len(compounds)} compounds")
        print(f"   First compound: {compounds[0].name}")
        print(f"   SMILES length: {len(compounds[0].smiles)}")
    except Exception as e:
        print(f"❌ Compound loading error: {e}")
        return False
    
    # Test 3: DeepDTA Adapter
    print("\n🧠 Test 3: DeepDTA Adapter")
    try:
        deepdta = DeepDTAAdapter()
        await deepdta.setup()
        
        # Test prediction
        test_protein = 'MKWVTFISLLFLFSSAYSGVFRRDTHKSEIAHRFKDLGEQFKQHLMNVKTRGQHLLGSLSTCPAGEDLLGKVNEFYANAHQSLVQASQPDGFKVLMSADNFMADLSEGTCVMQSGIRKYLNHMEGDYPELCDLAQAFEDLTHLDTEYSEFGT'
        result = await deepdta.predict_affinity(
            smiles=compounds[0].smiles,
            protein_sequence=test_protein,
            compound_id=compounds[0].compound_id,
            target_id='test_protein'
        )
        
        print(f"✅ DeepDTA prediction successful")
        print(f"   Predicted affinity: {result.predicted_affinity:.3f}")
        print(f"   Model version: {result.model_version}")
    except Exception as e:
        print(f"❌ DeepDTA error: {e}")
        return False
    
    # Test 4: Vina Binary
    print("\n⚗️ Test 4: Vina Binary")
    try:
        import subprocess
        import shutil
        
        # Check if vina_custom is in PATH
        vina_path = shutil.which('vina_custom')
        if vina_path:
            print(f"✅ Vina binary found: {vina_path}")
            
            # Test version (with library path)
            env = os.environ.copy()
            env['DYLD_LIBRARY_PATH'] = f"{os.environ.get('CONDA_PREFIX', '')}/lib:{env.get('DYLD_LIBRARY_PATH', '')}"
            
            result = subprocess.run(['vina_custom', '--version'], 
                                  capture_output=True, text=True, env=env)
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
    
    # Test 5: Environment
    print("\n🐍 Test 5: Environment")
    try:
        import numpy as np
        import pandas as pd
        import pydantic
        import httpx
        
        print(f"✅ Python: {sys.version.split()[0]}")
        print(f"✅ NumPy: {np.__version__}")
        print(f"✅ Pandas: {pd.__version__}")
        print(f"✅ Pydantic: {pydantic.__version__}")
        print(f"✅ HTTPX: {httpx.__version__}")
        
        # Check conda environment
        conda_prefix = os.environ.get('CONDA_PREFIX')
        if conda_prefix:
            print(f"✅ Conda environment: {Path(conda_prefix).name}")
        else:
            print("⚠️ No conda environment detected")
    except Exception as e:
        print(f"❌ Environment error: {e}")
        return False
    
    # Summary
    print("\n📊 Test Summary")
    print("=" * 50)
    print("✅ Configuration loading: WORKING")
    print("✅ Compound loading: WORKING")
    print("✅ DeepDTA adapter: WORKING")
    print("✅ Vina binary: WORKING")
    print("✅ Environment: WORKING")
    print("\n🎯 Status: READY FOR MCP SERVER INTEGRATION")
    print("   Next step: Start MCP servers with 'docker compose up -d'")
    
    return True

if __name__ == "__main__":
    success = asyncio.run(test_current_status())
    sys.exit(0 if success else 1)
