#!/usr/bin/env python3
"""
Interactive demo showing how to test each functionality step by step.

This script demonstrates practical testing approaches for the drug discovery pipeline.
"""

import asyncio
import json
import time
from pathlib import Path
import sys

# Add orchestrator to path
sys.path.append(str(Path(__file__).parent.parent))


class PipelineDemo:
    """Interactive demo of pipeline testing."""
    
    def __init__(self):
        self.step_counter = 1
    
    def step(self, title: str):
        """Print step header."""
        print(f"\n{'='*60}")
        print(f"STEP {self.step_counter}: {title}")
        print(f"{'='*60}")
        self.step_counter += 1
    
    def wait_for_user(self, message: str = "Press Enter to continue..."):
        """Wait for user input."""
        input(f"\n👆 {message}")
    
    async def demo_1_environment_check(self):
        """Demo 1: Check environment and dependencies."""
        self.step("Environment and Dependencies Check")
        
        print("🔍 Checking Python environment...")
        print(f"Python version: {sys.version}")
        
        try:
            import orchestrator
            print("✅ Orchestrator package importable")
        except ImportError as e:
            print(f"❌ Cannot import orchestrator: {e}")
            return False
        
        print("\n🔍 Checking required imports...")
        required_imports = [
            ("orchestrator.settings", "settings"),
            ("orchestrator.schemas.config", "RunConfig"),
            ("orchestrator.pipeline", "DrugDiscoveryPipeline"),
            ("orchestrator.mcp_clients", "KEGGClient"),
            ("orchestrator.adapters", "DeepDTAAdapter"),
            ("orchestrator.scoring", "ScoreNormalizer")
        ]
        
        for module, item in required_imports:
            try:
                exec(f"from {module} import {item}")
                print(f"✅ {module}.{item}")
            except ImportError as e:
                print(f"❌ {module}.{item}: {e}")
        
        self.wait_for_user()
        return True
    
    async def demo_2_create_test_data(self):
        """Demo 2: Create test data for demonstrations."""
        self.step("Creating Test Data")
        
        print("📁 Creating test compound files...")
        
        # Create compounds directory
        compounds_dir = Path("data/compounds")
        compounds_dir.mkdir(exist_ok=True, parents=True)
        
        # Create small test set
        small_test = compounds_dir / "demo_small.smi"
        with open(small_test, "w") as f:
            compounds = [
                "CCO ethanol",
                "CC(C)O isopropanol", 
                "C1CCCCC1 cyclohexane",
                "CC(=O)O acetic_acid",
                "C1=CC=CC=C1 benzene"
            ]
            for compound in compounds:
                f.write(f"{compound}\n")
        
        print(f"✅ Created {small_test} with {len(compounds)} compounds")
        
        # Create larger test set
        medium_test = compounds_dir / "demo_medium.smi"
        with open(medium_test, "w") as f:
            # Add the small set plus more
            for compound in compounds:
                f.write(f"{compound}\n")
            
            additional = [
                "CCN(CC)CC diethylamine",
                "CC(C)(C)O tert-butanol",
                "C1CCOC1 tetrahydrofuran",
                "CCO[H] ethanol_h",
                "C1=CC=C(C=C1)O phenol",
                "CC(C)CC(C)C isohexane",
                "CCCCO 1-butanol",
                "CC1=CC=CC=C1 toluene",
                "ClCCCl 1,3-dichloropropane",
                "CCCCCC hexane"
            ]
            for compound in additional:
                f.write(f"{compound}\n")
        
        print(f"✅ Created {medium_test} with {len(compounds + additional)} compounds")
        
        # Create test configurations
        configs_dir = Path("configs")
        configs_dir.mkdir(exist_ok=True, parents=True)
        
        # Small test config
        small_config = {
            "disease_terms": ["lung cancer"],
            "organism": "Homo sapiens",
            "compounds": {
                "input_paths": ["data/compounds/demo_small.smi"],
                "max_compounds": 5
            },
            "max_targets": 3,
            "scoring": {
                "weights": {
                    "deepdta": 0.6,
                    "docking": 0.3,
                    "evidence": 0.1
                }
            },
            "output_dir": "data/outputs/demo",
            "debug_mode": True
        }
        
        import yaml
        with open(configs_dir / "demo_small.yaml", "w") as f:
            yaml.dump(small_config, f, default_flow_style=False)
        
        print("✅ Created demo_small.yaml configuration")
        
        self.wait_for_user()
    
    async def demo_3_test_individual_components(self):
        """Demo 3: Test individual components."""
        self.step("Testing Individual Components")
        
        print("🔧 Testing score normalization...")
        from orchestrator.scoring import ScoreNormalizer
        
        normalizer = ScoreNormalizer()
        test_scores = [1.0, 2.5, 3.0, 4.5, 5.0]
        
        print(f"Original scores: {test_scores}")
        z_normalized = normalizer.z_score_normalize(test_scores)
        print(f"Z-normalized: {[f'{x:.2f}' for x in z_normalized]}")
        
        self.wait_for_user("Continue to test MCP clients...")
        
        print("\n🔌 Testing MCP client initialization...")
        from orchestrator.mcp_clients import KEGGClient, UniProtClient
        
        # Test KEGG client
        print("Testing KEGG client...")
        kegg = KEGGClient()
        print(f"  URL: {kegg.base_url}")
        print(f"  Timeout: {kegg.timeout}s")
        
        try:
            connected = await kegg.test_connection()
            print(f"  Connection: {'✅ Success' if connected else '❌ Failed'}")
            
            if connected:
                print("  Testing disease pathway search...")
                pathways = await kegg.search_disease_pathways(["cancer"], max_pathways_per_disease=3)
                print(f"  Found {len(pathways)} pathways")
                
                if pathways:
                    pathway = pathways[0]
                    print(f"  Example: {pathway.get('pathway_id')} - {pathway.get('pathway_name', 'N/A')[:40]}...")
                    
        except Exception as e:
            print(f"  ⚠️ KEGG test failed (expected if service not running): {e}")
        
        finally:
            await kegg.close()
        
        self.wait_for_user("Continue to test adapters...")
        
        print("\n🤖 Testing ML adapter initialization...")
        from orchestrator.adapters import DeepDTAAdapter, GeminiMolAdapter
        
        # Test DeepDTA adapter
        print("Testing DeepDTA adapter...")
        deepdta = DeepDTAAdapter()
        print(f"  Model directory: {deepdta.deepdta_dir}")
        
        try:
            setup_ok = await deepdta.setup()
            print(f"  Setup: {'✅ Success' if setup_ok else '❌ Failed'}")
            
            if setup_ok:
                # Test encoding
                smiles = "CCO"
                sequence = "MVLSPADKTNVKAAW"
                
                compound_encoded = deepdta._encode_compound(smiles)
                protein_encoded = deepdta._encode_protein(sequence)
                
                print(f"  Encoded compound length: {len(compound_encoded)}")
                print(f"  Encoded protein length: {len(protein_encoded)}")
                
                # Test prediction
                score = await deepdta.predict_affinity(smiles, sequence)
                print(f"  Predicted affinity: {score.predicted_affinity:.3f}")
                
        except Exception as e:
            print(f"  ⚠️ DeepDTA test failed (expected if models not available): {e}")
        
        self.wait_for_user()
    
    async def demo_4_test_pipeline_stages(self):
        """Demo 4: Test pipeline stages individually."""
        self.step("Testing Pipeline Stages")
        
        print("🔄 Testing pipeline components...")
        from orchestrator.pipeline import DrugDiscoveryPipeline
        from orchestrator.schemas.config import RunConfig
        
        pipeline = DrugDiscoveryPipeline()
        
        # Load configuration
        config = RunConfig(
            disease_terms=["lung cancer"],
            compounds={"input_paths": ["data/compounds/demo_small.smi"]},
            max_targets=2
        )
        
        print("Configuration loaded:")
        print(f"  Disease terms: {config.disease_terms}")
        print(f"  Compound files: {config.compounds.input_paths}")
        print(f"  Max targets: {config.max_targets}")
        
        self.wait_for_user("Test compound loading...")
        
        # Test compound loading
        print("\n📦 Testing compound loading...")
        
        class MockRun:
            def add_stage_result(self, result):
                print(f"    Stage result: {result.status.value}")
        
        mock_run = MockRun()
        
        try:
            compounds = await pipeline._load_compounds(config, mock_run)
            print(f"✅ Loaded {len(compounds)} compounds")
            
            for i, compound in enumerate(compounds):
                print(f"  {i+1}. {compound.compound_id}: {compound.smiles} ({compound.name})")
                
        except Exception as e:
            print(f"❌ Compound loading failed: {e}")
        
        self.wait_for_user("Test target discovery...")
        
        # Test target discovery
        print("\n🎯 Testing target discovery...")
        try:
            targets = await pipeline._discover_targets(config, mock_run)
            print(f"✅ Discovered {len(targets)} targets")
            
            for i, target in enumerate(targets[:3]):  # Show first 3
                print(f"  {i+1}. {target.target_id}: {target.gene_name}")
                print(f"      Sources: {target.discovery_source}")
                
        except Exception as e:
            print(f"⚠️ Target discovery failed (expected if MCP services not running): {e}")
        
        self.wait_for_user()
    
    async def demo_5_test_scoring_integration(self):
        """Demo 5: Test scoring system integration."""
        self.step("Testing Scoring System Integration")
        
        print("📊 Creating mock scoring data...")
        
        # Create realistic mock data
        mock_pairs = [
            {
                "compound_id": "comp_001",
                "target_id": "EGFR",
                "compound_smiles": "CCO",
                "target_gene": "EGFR",
                "deepdta_score": 8.5,
                "docking_score": -7.2,
                "evidence_score": 6.0
            },
            {
                "compound_id": "comp_002", 
                "target_id": "EGFR",
                "compound_smiles": "CC(C)O",
                "target_gene": "EGFR",
                "deepdta_score": 7.8,
                "docking_score": -8.1,
                "evidence_score": 5.5
            },
            {
                "compound_id": "comp_001",
                "target_id": "TP53",
                "compound_smiles": "CCO",
                "target_gene": "TP53",
                "deepdta_score": 6.9,
                "docking_score": -6.8,
                "evidence_score": 7.2
            },
            {
                "compound_id": "comp_003",
                "target_id": "KRAS",
                "compound_smiles": "C1CCCCC1",
                "target_gene": "KRAS",
                "deepdta_score": 9.1,
                "docking_score": -7.5,
                "evidence_score": 4.8
            }
        ]
        
        print(f"Created {len(mock_pairs)} compound-target pairs")
        
        # Test score fusion
        print("\n🔧 Testing score fusion...")
        from orchestrator.scoring import ScoreFusion
        
        weights = {"deepdta": 0.6, "docking": 0.3, "evidence": 0.1}
        fusion = ScoreFusion(weights)
        
        print(f"Using weights: {weights}")
        
        combined_scores = fusion.combine_scores(mock_pairs, normalization_method="zscore")
        
        print("\n📈 Ranked results:")
        for i, result in enumerate(combined_scores[:3]):  # Top 3
            print(f"  Rank {result['rank']}: {result['compound_id']}-{result['target_id']}")
            print(f"    Combined score: {result['combined_score']:.3f}")
            print(f"    DeepDTA: {result['deepdta_score']} → {result['deepdta_normalized']:.3f}")
            print(f"    Docking: {result['docking_score']} → {result['docking_normalized']:.3f}")
            print(f"    Evidence: {result['evidence_score']} → {result['evidence_normalized']:.3f}")
            print()
        
        self.wait_for_user()
    
    async def demo_6_test_configuration_validation(self):
        """Demo 6: Test configuration validation."""
        self.step("Testing Configuration Validation")
        
        print("⚙️ Testing configuration validation...")
        from orchestrator.schemas.config import RunConfig
        
        # Test valid configuration
        print("\n✅ Testing valid configuration...")
        try:
            valid_config = RunConfig(
                disease_terms=["lung cancer"],
                compounds={
                    "input_paths": ["data/compounds/demo_small.smi"]
                },
                scoring={
                    "weights": {
                        "deepdta": 0.6,
                        "docking": 0.3,
                        "evidence": 0.1
                    }
                }
            )
            print("Valid configuration created successfully")
            print(f"  Disease terms: {valid_config.disease_terms}")
            print(f"  Organism: {valid_config.organism}")
            print(f"  Scoring weights: {valid_config.scoring.weights}")
            
        except Exception as e:
            print(f"❌ Valid configuration failed: {e}")
        
        self.wait_for_user("Test invalid configurations...")
        
        # Test invalid configurations
        print("\n❌ Testing invalid configurations...")
        
        invalid_configs = [
            {
                "name": "Empty disease terms",
                "config": {
                    "disease_terms": [],  # Invalid
                    "compounds": {"input_paths": ["test.smi"]}
                }
            },
            {
                "name": "Invalid scoring weights",
                "config": {
                    "disease_terms": ["cancer"],
                    "compounds": {"input_paths": ["test.smi"]},
                    "scoring": {
                        "weights": {
                            "deepdta": 0.5,
                            "docking": 0.3,
                            "evidence": 0.1  # Sum = 0.9, not 1.0
                        }
                    }
                }
            }
        ]
        
        for test_case in invalid_configs:
            print(f"\nTesting: {test_case['name']}")
            try:
                RunConfig(**test_case['config'])
                print(f"  ⚠️ Configuration was accepted (should have failed)")
            except Exception as e:
                print(f"  ✅ Configuration properly rejected: {str(e)[:80]}...")
        
        self.wait_for_user()
    
    async def demo_7_end_to_end_dry_run(self):
        """Demo 7: End-to-end dry run test."""
        self.step("End-to-End Dry Run Test")
        
        print("🏃‍♂️ Testing complete pipeline in dry run mode...")
        
        # Use CLI interface for dry run
        print("\n🖥️ Using CLI interface...")
        import subprocess
        import os
        
        # Change to project directory
        original_dir = os.getcwd()
        project_dir = Path(__file__).parent.parent
        os.chdir(project_dir)
        
        try:
            # Test configuration validation
            print("Testing configuration validation...")
            result = subprocess.run([
                sys.executable, "-m", "orchestrator.cli", 
                "validate", "configs/demo_small.yaml"
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                print("✅ Configuration validation passed")
                print(f"Output: {result.stdout}")
            else:
                print("❌ Configuration validation failed")
                print(f"Error: {result.stderr}")
            
            self.wait_for_user("Run dry run test...")
            
            # Create dry run config
            dry_run_config = Path("configs/demo_small.yaml")
            if dry_run_config.exists():
                import yaml
                with open(dry_run_config) as f:
                    config = yaml.safe_load(f)
                config["dry_run"] = True
                
                with open("configs/demo_dry_run.yaml", "w") as f:
                    yaml.dump(config, f, default_flow_style=False)
                
                print("Running dry run...")
                result = subprocess.run([
                    sys.executable, "-m", "orchestrator.cli",
                    "run", "configs/demo_dry_run.yaml"
                ], capture_output=True, text=True, timeout=30)
                
                if result.returncode == 0:
                    print("✅ Dry run completed successfully")
                    print("Output snippet:", result.stdout[-200:] if result.stdout else "No output")
                else:
                    print("❌ Dry run failed")
                    print("Error:", result.stderr[-200:] if result.stderr else "No error output")
        
        except subprocess.TimeoutExpired:
            print("⚠️ Dry run timed out (expected for complex operations)")
        except Exception as e:
            print(f"❌ Dry run test failed: {e}")
        finally:
            os.chdir(original_dir)
        
        self.wait_for_user()
    
    def demo_8_summary_and_next_steps(self):
        """Demo 8: Summary and next steps."""
        self.step("Summary and Next Steps")
        
        print("🎯 Testing Summary")
        print("-" * 30)
        print("✅ Environment check completed")
        print("✅ Test data created")
        print("✅ Individual components tested")
        print("✅ Pipeline stages tested")
        print("✅ Scoring integration tested")
        print("✅ Configuration validation tested")
        print("✅ End-to-end dry run tested")
        
        print("\n🚀 Next Steps for Full Testing")
        print("-" * 35)
        print("1. Start MCP services:")
        print("   docker compose up -d")
        print()
        print("2. Run comprehensive tests:")
        print("   python scripts/test_functionality.py")
        print()
        print("3. Test individual components:")
        print("   python scripts/test_component.py mcp kegg")
        print("   python scripts/test_component.py adapter deepdta")
        print()
        print("4. Run quick functionality check:")
        print("   ./scripts/quick_test.sh")
        print()
        print("5. Try a real pipeline run:")
        print("   python -m orchestrator.cli run configs/demo_small.yaml")
        print()
        print("6. Use the API:")
        print("   curl -X POST http://localhost:8000/runs -d @config.json")
        print()
        print("📚 For detailed testing documentation, see: TESTING.md")


async def main():
    """Run the interactive demo."""
    demo = PipelineDemo()
    
    print("🧪 MCP Drug Discovery Pipeline - Interactive Testing Demo")
    print("=" * 60)
    print("This demo will walk you through testing each component of the pipeline.")
    print("You can stop at any point with Ctrl+C.")
    
    try:
        # Run demo steps
        await demo.demo_1_environment_check()
        await demo.demo_2_create_test_data()
        await demo.demo_3_test_individual_components()
        await demo.demo_4_test_pipeline_stages()
        await demo.demo_5_test_scoring_integration()
        await demo.demo_6_test_configuration_validation()
        await demo.demo_7_end_to_end_dry_run()
        demo.demo_8_summary_and_next_steps()
        
        print("\n🎉 Demo completed successfully!")
        print("You now know how to test each component of the pipeline.")
        
    except KeyboardInterrupt:
        print("\n👋 Demo interrupted by user. Goodbye!")
    except Exception as e:
        print(f"\n💥 Demo failed with error: {e}")
        raise


if __name__ == "__main__":
    asyncio.run(main())
