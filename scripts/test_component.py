#!/usr/bin/env python3
"""
Individual component testing script.

Usage:
    python scripts/test_component.py mcp kegg
    python scripts/test_component.py adapter geminimol
    python scripts/test_component.py pipeline compound_loading
    python scripts/test_component.py scoring normalization
"""

import asyncio
import sys
import argparse
from pathlib import Path

# Add orchestrator to path
sys.path.append(str(Path(__file__).parent.parent))


async def test_mcp_client(client_name: str):
    """Test specific MCP client."""
    print(f"🔌 Testing {client_name.upper()} MCP Client")
    print("=" * 40)
    
    try:
        if client_name.lower() == "kegg":
            from orchestrator.mcp_clients import KEGGClient
            client = KEGGClient()
            
            print(f"URL: {client.base_url}")
            connected = await client.test_connection()
            print(f"Connection: {'✅ Connected' if connected else '❌ Failed'}")
            
            if connected:
                pathways = await client.search_disease_pathways(["lung cancer"])
                print(f"Disease pathway search: Found {len(pathways)} pathways")
                
                if pathways:
                    pathway = pathways[0]
                    print(f"Example pathway: {pathway.get('pathway_id')} - {pathway.get('pathway_name', 'N/A')[:50]}")
            
            await client.close()
            
        elif client_name.lower() == "uniprot":
            from orchestrator.mcp_clients import UniProtClient
            client = UniProtClient()
            
            print(f"URL: {client.base_url}")
            connected = await client.test_connection()
            print(f"Connection: {'✅ Connected' if connected else '❌ Failed'}")
            
            if connected:
                proteins = await client.search_proteins("insulin", limit=3)
                print(f"Protein search: Found {len(proteins)} proteins")
                
                if proteins:
                    protein = proteins[0]
                    print(f"Example protein: {protein.get('accession')} - {protein.get('protein_name', 'N/A')[:50]}")
            
            await client.close()
            
        elif client_name.lower() == "string":
            from orchestrator.mcp_clients import STRINGClient
            client = STRINGClient()
            
            print(f"URL: {client.base_url}")
            connected = await client.test_connection()
            print(f"Connection: {'✅ Connected' if connected else '❌ Failed'}")
            
            if connected:
                interactions = await client.get_protein_interactions(["ENSP00000000233"])
                print(f"PPI search: Found {len(interactions)} interactions")
            
            await client.close()
            
        else:
            print(f"❌ Unknown MCP client: {client_name}")
            print("Available clients: kegg, uniprot, string, reactome, proteinatlas, pdb, chembl")
            
    except Exception as e:
        print(f"❌ Error testing {client_name}: {e}")


async def test_adapter(adapter_name: str):
    """Test specific ML adapter."""
    print(f"🤖 Testing {adapter_name.upper()} Adapter")
    print("=" * 40)
    
    try:
        if adapter_name.lower() == "geminimol":
            from orchestrator.adapters import GeminiMolAdapter
            adapter = GeminiMolAdapter()
            
            print(f"GeminiMol directory: {adapter.geminimol_dir}")
            
            setup_ok = await adapter.setup()
            print(f"Setup: {'✅ Success' if setup_ok else '❌ Failed'}")
            
            if setup_ok:
                print("Testing embedding computation...")
                embedding = await adapter.compute_embedding("CCO")
                print(f"Embedding size: {len(embedding)}")
                print(f"Sample values: {embedding[:5]}")
                
                # Test batch computation
                smiles_list = ["CCO", "CC(C)O", "C1CCCCC1"]
                embeddings = await adapter.compute_batch_embeddings(smiles_list)
                print(f"Batch embeddings: {len(embeddings)} compounds processed")
            
        elif adapter_name.lower() == "vina":
            from orchestrator.adapters import VinaAdapter
            adapter = VinaAdapter()
            
            print(f"Vina binary: {adapter.vina_binary}")
            print(f"Exhaustiveness: {adapter.exhaustiveness}")
            print(f"Number of modes: {adapter.num_modes}")
            
            setup_ok = await adapter.setup()
            print(f"Setup: {'✅ Success' if setup_ok else '❌ Failed'}")
            
            # Test ligand preparation
            print("Testing ligand preparation...")
            ligand_path = await adapter._prepare_ligand("CCO")
            if ligand_path:
                print(f"Ligand file: {ligand_path}")
                print(f"File exists: {Path(ligand_path).exists()}")
            
        elif adapter_name.lower() == "ouroboros":
            from orchestrator.adapters import OuroborosAdapter
            adapter = OuroborosAdapter()
            
            setup_ok = await adapter.setup()
            print(f"Setup: {'✅ Success' if setup_ok else '❌ Failed'}")
            
            if setup_ok:
                print("Testing feature unification...")
                geminimol_features = [0.1, 0.2, 0.3, 0.4, 0.5]
                additional_features = {"rdkit": [0.6, 0.7], "other": [0.8]}
                
                unified = await adapter.unify_features(geminimol_features, additional_features)
                print(f"Original size: {len(geminimol_features)}")
                print(f"Unified size: {len(unified)}")
                
                # Test similarity
                similarity = await adapter.compute_similarity(unified, unified)
                print(f"Self-similarity: {similarity:.3f} (should be ~1.0)")
        
        else:
            print(f"❌ Unknown adapter: {adapter_name}")
            print("Available adapters: geminimol, vina, ouroboros")
            
    except Exception as e:
        print(f"❌ Error testing {adapter_name}: {e}")


async def test_pipeline_stage(stage_name: str):
    """Test specific pipeline stage."""
    print(f"🔄 Testing Pipeline Stage: {stage_name.upper()}")
    print("=" * 40)
    
    try:
        from orchestrator.pipeline import DrugDiscoveryPipeline
        from orchestrator.schemas.config import RunConfig
        
        pipeline = DrugDiscoveryPipeline()
        
        # Create test data
        test_compounds_file = Path("data/compounds/test_component.smi")
        test_compounds_file.parent.mkdir(exist_ok=True, parents=True)
        
        with open(test_compounds_file, "w") as f:
            f.write("CCO ethanol\n")
            f.write("CC(C)O isopropanol\n")
            f.write("C1CCCCC1 cyclohexane\n")
        
        config = RunConfig(
            disease_terms=["lung cancer"],
            compounds={"input_paths": [str(test_compounds_file)]},
            max_targets=3
        )
        
        class MockRun:
            def add_stage_result(self, result): 
                print(f"Stage result: {result.status.value}")
        
        mock_run = MockRun()
        
        if stage_name.lower() == "compound_loading":
            print("Testing compound loading...")
            compounds = await pipeline._load_compounds(config, mock_run)
            print(f"Loaded compounds: {len(compounds)}")
            
            for i, compound in enumerate(compounds):
                print(f"  {i+1}. {compound.compound_id}: {compound.smiles} ({compound.name})")
        
        elif stage_name.lower() == "target_discovery":
            print("Testing target discovery...")
            targets = await pipeline._discover_targets(config, mock_run)
            print(f"Discovered targets: {len(targets)}")
            
            for i, target in enumerate(targets[:5]):  # Show first 5
                print(f"  {i+1}. {target.target_id}: {target.gene_name}")
                
        else:
            print(f"❌ Unknown pipeline stage: {stage_name}")
            print("Available stages: compound_loading, target_discovery")
            
    except Exception as e:
        print(f"❌ Error testing pipeline stage {stage_name}: {e}")


def test_scoring(test_type: str):
    """Test scoring system components."""
    print(f"📊 Testing Scoring: {test_type.upper()}")
    print("=" * 40)
    
    try:
        if test_type.lower() == "normalization":
            from orchestrator.scoring import ScoreNormalizer
            
            normalizer = ScoreNormalizer()
            test_scores = [1.0, 2.5, 3.0, 4.5, 5.0, 2.1, 3.7]
            
            print(f"Original scores: {test_scores}")
            
            # Test z-score normalization
            z_normalized = normalizer.z_score_normalize(test_scores)
            print(f"Z-score normalized: {[f'{x:.3f}' for x in z_normalized]}")
            
            # Test min-max normalization
            minmax_normalized = normalizer.min_max_normalize(test_scores)
            print(f"Min-max normalized: {[f'{x:.3f}' for x in minmax_normalized]}")
            
            # Test rank normalization
            rank_normalized = normalizer.rank_normalize(test_scores)
            print(f"Rank normalized: {[f'{x:.3f}' for x in rank_normalized]}")
            
        elif test_type.lower() == "fusion":
            from orchestrator.scoring import ScoreFusion
            
            weights = {"similarity": 0.5, "pharmacophore": 0.2, "docking": 0.1, "evidence": 0.2}
            fusion = ScoreFusion(weights)
            
            score_data = [
                {"pair_id": "comp1-targ1", "docking_score": -7.2, "evidence_score": 6.0},
                {"pair_id": "comp1-targ2", "docking_score": -8.1, "evidence_score": 5.5},
                {"pair_id": "comp2-targ1", "docking_score": -6.8, "evidence_score": 7.2},
                {"pair_id": "comp2-targ2", "docking_score": -7.5, "evidence_score": 4.8}
            ]
            
            print("Original score data:")
            for item in score_data:
                print(f"  {item['pair_id']}: Docking={item['docking_score']}, Evidence={item['evidence_score']}")
            
            combined = fusion.combine_scores(score_data)
            
            print("\nCombined and ranked scores:")
            for item in combined:
                print(f"  Rank {item['rank']}: {item['pair_id']} = {item['combined_score']:.3f}")
                
        else:
            print(f"❌ Unknown scoring test: {test_type}")
            print("Available scoring tests: normalization, fusion")
            
    except Exception as e:
        print(f"❌ Error testing scoring {test_type}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Test individual components of the MCP Drug Discovery Pipeline")
    parser.add_argument("component", choices=["mcp", "adapter", "pipeline", "scoring"], 
                       help="Component type to test")
    parser.add_argument("name", help="Specific component name to test")
    
    args = parser.parse_args()
    
    print(f"🧪 Component Testing: {args.component.upper()} - {args.name.upper()}")
    print("=" * 60)
    
    try:
        if args.component == "mcp":
            asyncio.run(test_mcp_client(args.name))
        elif args.component == "adapter":
            asyncio.run(test_adapter(args.name))
        elif args.component == "pipeline":
            asyncio.run(test_pipeline_stage(args.name))
        elif args.component == "scoring":
            test_scoring(args.name)
        
        print("\n✅ Component test completed!")
        
    except KeyboardInterrupt:
        print("\n❌ Test interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Test failed with error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
