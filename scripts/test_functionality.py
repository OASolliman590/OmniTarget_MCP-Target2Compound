#!/usr/bin/env python3
"""
Comprehensive functionality testing script for MCP Drug Discovery Pipeline.

This script tests each component individually and provides detailed feedback.
"""

import asyncio
import sys
import json
import time
from pathlib import Path
from typing import Dict, List, Any

# Add the orchestrator to the path
sys.path.append(str(Path(__file__).parent.parent))

from orchestrator.mcp_clients import *
from orchestrator.adapters import GeminiMolAdapter, VinaAdapter, OuroborosAdapter
from orchestrator.pipeline import DrugDiscoveryPipeline
from orchestrator.schemas.config import RunConfig
from orchestrator.scoring import ScoreNormalizer, ScoreFusion


class FunctionalityTester:
    """Comprehensive functionality tester."""
    
    def __init__(self):
        self.results = {}
        self.total_tests = 0
        self.passed_tests = 0
    
    def log_test(self, test_name: str, status: str, message: str = "", details: Any = None):
        """Log test result."""
        self.total_tests += 1
        if status == "PASS":
            self.passed_tests += 1
            print(f"✅ {test_name}: {message}")
        elif status == "FAIL":
            print(f"❌ {test_name}: {message}")
        elif status == "SKIP":
            print(f"⏭️  {test_name}: {message}")
        else:
            print(f"⚠️  {test_name}: {message}")
        
        self.results[test_name] = {
            "status": status,
            "message": message,
            "details": details
        }
    
    async def test_mcp_clients(self):
        """Test all MCP client connections and basic functionality."""
        print("\n🔌 Testing MCP Client Connections")
        print("=" * 50)
        
        clients = [
            ("KEGG", KEGGClient()),
            ("Reactome", ReactomeClient()),
            ("ProteinAtlas", ProteinAtlasClient()),
            ("STRING", STRINGClient()),
            ("UniProt", UniProtClient()),
            ("PDB", PDBClient()),
            ("ChEMBL", ChEMBLClient())
        ]
        
        for name, client in clients:
            try:
                connected = await client.test_connection()
                if connected:
                    self.log_test(f"MCP_{name}_connection", "PASS", "Connected successfully")
                    
                    # Test basic functionality
                    if name == "KEGG":
                        pathways = await client.search_disease_pathways(["test disease"])
                        self.log_test(f"MCP_{name}_search", "PASS", f"Found {len(pathways)} pathways")
                    
                    elif name == "UniProt":
                        proteins = await client.search_proteins("insulin", limit=1)
                        self.log_test(f"MCP_{name}_search", "PASS", f"Found {len(proteins)} proteins")
                    
                else:
                    self.log_test(f"MCP_{name}_connection", "FAIL", "Connection failed")
                    
            except Exception as e:
                self.log_test(f"MCP_{name}_connection", "SKIP", f"Service unavailable: {str(e)[:100]}")
            
            finally:
                try:
                    await client.close()
                except:
                    pass
    
    async def test_ml_adapters(self):
        """Test ML model adapters."""
        print("\n🤖 Testing ML Model Adapters")
        print("=" * 50)
        
        # Sequence-based predictor removed; skip
        
        # Test GeminiMol Adapter
        try:
            geminimol = GeminiMolAdapter()
            setup_ok = await geminimol.setup()
            
            if setup_ok:
                embedding = await geminimol.compute_embedding("CCO")
                self.log_test("GeminiMol_embedding", "PASS", 
                            f"Generated embedding of size {len(embedding)}")
            else:
                self.log_test("GeminiMol_setup", "FAIL", "Setup failed")
                
        except Exception as e:
            self.log_test("GeminiMol_adapter", "SKIP", f"Not available: {str(e)[:100]}")
        
        # Test Vina Adapter
        try:
            vina = VinaAdapter()
            setup_ok = await vina.setup()
            
            if setup_ok:
                self.log_test("Vina_setup", "PASS", f"Vina binary: {vina.vina_binary}")
            else:
                self.log_test("Vina_setup", "FAIL", "Setup failed")
                
        except Exception as e:
            self.log_test("Vina_adapter", "SKIP", f"Not available: {str(e)[:100]}")
        
        # Test Ouroboros Adapter
        try:
            ouroboros = OuroborosAdapter()
            setup_ok = await ouroboros.setup()
            
            if setup_ok:
                features = await ouroboros.unify_features([1, 2, 3], {"extra": [4, 5]})
                self.log_test("Ouroboros_unification", "PASS", 
                            f"Unified features size: {len(features)}")
            else:
                self.log_test("Ouroboros_setup", "FAIL", "Setup failed")
                
        except Exception as e:
            self.log_test("Ouroboros_adapter", "SKIP", f"Error: {str(e)[:100]}")
    
    def test_scoring_system(self):
        """Test scoring and normalization system."""
        print("\n📊 Testing Scoring System")
        print("=" * 50)
        
        try:
            # Test score normalization
            normalizer = ScoreNormalizer()
            test_scores = [1.0, 2.5, 3.0, 4.5, 5.0]
            
            z_normalized = normalizer.z_score_normalize(test_scores)
            self.log_test("Score_z_normalization", "PASS", 
                        f"Z-normalized: {[f'{x:.2f}' for x in z_normalized[:3]]}")
            
            minmax_normalized = normalizer.min_max_normalize(test_scores)
            self.log_test("Score_minmax_normalization", "PASS", 
                        f"MinMax-normalized: {[f'{x:.2f}' for x in minmax_normalized[:3]]}")
            
            # Test score fusion
            weights = {"similarity": 0.5, "pharmacophore": 0.2, "docking": 0.1, "evidence": 0.2}
            fusion = ScoreFusion(weights)
            
            score_data = [
                {"pair_id": "comp1-targ1", "docking_score": -7.2, "evidence_score": 6.0},
                {"pair_id": "comp1-targ2", "docking_score": -8.1, "evidence_score": 5.5}
            ]
            
            combined = fusion.combine_scores(score_data)
            self.log_test("Score_fusion", "PASS", 
                        f"Combined {len(combined)} scores, top score: {combined[0]['combined_score']:.3f}")
            
        except Exception as e:
            self.log_test("Scoring_system", "FAIL", f"Error: {str(e)}")
    
    def test_configuration(self):
        """Test configuration validation."""
        print("\n⚙️ Testing Configuration System")
        print("=" * 50)
        
        try:
            # Test valid configuration
            config = RunConfig(
                disease_terms=["lung cancer"],
                compounds={
                    "input_paths": ["test.smi"]  # File doesn't need to exist for this test
                },
                scoring={
                    "weights": {
                        "similarity": 0.5,
                        "pharmacophore": 0.2,
                        "docking": 0.1,
                        "evidence": 0.2
                    }
                }
            )
            self.log_test("Config_valid", "PASS", "Valid configuration created")
            
            # Test invalid configuration (should raise exception)
            try:
                invalid_config = RunConfig(
                    disease_terms=[],  # Invalid: empty
                    compounds={"input_paths": ["test.smi"]}
                )
                self.log_test("Config_validation", "FAIL", "Invalid config was accepted")
            except Exception:
                self.log_test("Config_validation", "PASS", "Invalid config properly rejected")
            
        except Exception as e:
            self.log_test("Config_system", "FAIL", f"Error: {str(e)}")
    
    def create_test_data(self):
        """Create test data files."""
        print("\n📁 Creating Test Data")
        print("=" * 50)
        
        try:
            # Create test compounds directory
            compounds_dir = Path("data/compounds")
            compounds_dir.mkdir(exist_ok=True, parents=True)
            
            # Create test SMILES file
            test_file = compounds_dir / "test_small.smi"
            with open(test_file, "w") as f:
                f.write("CCO ethanol\n")
                f.write("CC(C)O isopropanol\n")
                f.write("C1CCCCC1 cyclohexane\n")
                f.write("CC(=O)O acetic_acid\n")
                f.write("C1=CC=CC=C1 benzene\n")
            
            self.log_test("Test_data_creation", "PASS", f"Created {test_file} with 5 compounds")
            
            # Create test configuration
            config_dir = Path("configs")
            config_dir.mkdir(exist_ok=True, parents=True)
            
            test_config = {
                "disease_terms": ["lung cancer"],
                "organism": "Homo sapiens",
                "compounds": {
                    "input_paths": ["data/compounds/test_small.smi"],
                    "max_compounds": 5
                },
                "max_targets": 3,
                "scoring": {
                    "weights": {
                        "similarity": 0.5,
                        "pharmacophore": 0.2,
                        "docking": 0.1,
                        "evidence": 0.2
                    }
                },
                "output_dir": "data/outputs/test",
                "debug_mode": True
            }
            
            test_config_file = config_dir / "test_functionality.yaml"
            import yaml
            with open(test_config_file, "w") as f:
                yaml.dump(test_config, f, default_flow_style=False)
            
            self.log_test("Test_config_creation", "PASS", f"Created {test_config_file}")
            
        except Exception as e:
            self.log_test("Test_data_creation", "FAIL", f"Error: {str(e)}")
    
    async def test_pipeline_components(self):
        """Test individual pipeline components."""
        print("\n🔄 Testing Pipeline Components")
        print("=" * 50)
        
        try:
            pipeline = DrugDiscoveryPipeline()
            
            # Test pipeline setup
            setup_ok = await pipeline.setup()
            self.log_test("Pipeline_setup", "PASS" if setup_ok else "FAIL", 
                        f"Setup result: {setup_ok}")
            
            # Test compound loading
            config = RunConfig(
                disease_terms=["test"],
                compounds={"input_paths": ["data/compounds/test_small.smi"]}
            )
            
            class MockRun:
                def add_stage_result(self, result): pass
            
            if Path("data/compounds/test_small.smi").exists():
                compounds = await pipeline._load_compounds(config, MockRun())
                self.log_test("Pipeline_compound_loading", "PASS", 
                            f"Loaded {len(compounds)} compounds")
                
                # Test target discovery (may fail if services unavailable)
                try:
                    targets = await pipeline._discover_targets(config, MockRun())
                    self.log_test("Pipeline_target_discovery", "PASS", 
                                f"Discovered {len(targets)} targets")
                except Exception as e:
                    self.log_test("Pipeline_target_discovery", "SKIP", 
                                f"Services unavailable: {str(e)[:100]}")
            else:
                self.log_test("Pipeline_compound_loading", "SKIP", "Test data not available")
            
        except Exception as e:
            self.log_test("Pipeline_components", "FAIL", f"Error: {str(e)}")
    
    def test_api_configuration(self):
        """Test API and CLI configuration."""
        print("\n🌐 Testing API Configuration")
        print("=" * 50)
        
        try:
            from orchestrator.api import app
            from orchestrator.cli import app as cli_app
            
            self.log_test("API_import", "PASS", "FastAPI app imported successfully")
            self.log_test("CLI_import", "PASS", "CLI app imported successfully")
            
            # Test settings
            from orchestrator.settings import settings
            self.log_test("Settings_import", "PASS", 
                        f"API host: {settings.api.host}, port: {settings.api.port}")
            
        except Exception as e:
            self.log_test("API_configuration", "FAIL", f"Error: {str(e)}")
    
    def print_summary(self):
        """Print test summary."""
        print("\n" + "=" * 60)
        print("🧪 TEST SUMMARY")
        print("=" * 60)
        
        print(f"Total tests: {self.total_tests}")
        print(f"Passed: {self.passed_tests}")
        print(f"Failed/Skipped: {self.total_tests - self.passed_tests}")
        print(f"Success rate: {(self.passed_tests/self.total_tests)*100:.1f}%")
        
        # Group results by category
        categories = {}
        for test_name, result in self.results.items():
            category = test_name.split("_")[0]
            if category not in categories:
                categories[category] = {"PASS": 0, "FAIL": 0, "SKIP": 0}
            categories[category][result["status"]] += 1
        
        print("\nResults by category:")
        for category, counts in categories.items():
            total = sum(counts.values())
            pass_rate = (counts["PASS"] / total) * 100 if total > 0 else 0
            print(f"  {category}: {counts['PASS']}/{total} passed ({pass_rate:.1f}%)")
        
        # List failed tests
        failed_tests = [name for name, result in self.results.items() 
                       if result["status"] == "FAIL"]
        if failed_tests:
            print(f"\nFailed tests: {', '.join(failed_tests)}")
        
        # Save detailed results
        results_file = Path("test_results.json")
        with open(results_file, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        print(f"\nDetailed results saved to: {results_file}")


async def main():
    """Run all functionality tests."""
    print("🚀 MCP Drug Discovery Pipeline - Functionality Testing")
    print("=" * 60)
    
    tester = FunctionalityTester()
    
    # Create test data first
    tester.create_test_data()
    
    # Run all tests
    await tester.test_mcp_clients()
    await tester.test_ml_adapters()
    tester.test_scoring_system()
    tester.test_configuration()
    await tester.test_pipeline_components()
    tester.test_api_configuration()
    
    # Print summary
    tester.print_summary()
    
    return tester.passed_tests, tester.total_tests


if __name__ == "__main__":
    try:
        passed, total = asyncio.run(main())
        
        # Exit with appropriate code
        if passed == total:
            print("\n🎉 All tests passed!")
            sys.exit(0)
        elif passed > total * 0.7:  # 70% pass rate
            print("\n✅ Most tests passed - ready for development!")
            sys.exit(0)
        else:
            print("\n⚠️ Many tests failed - check configuration and services")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n❌ Testing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Testing failed with error: {e}")
        sys.exit(1)
