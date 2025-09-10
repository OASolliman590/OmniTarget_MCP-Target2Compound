const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock ChEMBL MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'ChEMBL MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/status', (req, res) => {
  res.json({
    status: 'running',
    name: 'ChEMBL MCP Server',
    description: 'Mock ChEMBL database server'
  });
});

app.get('/genes', (req, res) => {
  res.json({
    genes: [
      { id: 'CHEMBL123', name: 'Mock Compound', description: 'Mock chemical compound' },
      { id: 'CHEMBL456', name: 'Mock Drug', description: 'Mock drug molecule' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: 'CHEMBL123', name: 'Mock Compound', activity: 'Active' }
    ]
  });
});

// Add the missing /activities endpoint
app.get('/activities', (req, res) => {
  const target = req.query.target;
  const compound = req.query.compound;
  
  // Mock bioactivity data for different targets
  const mockActivities = {
    'CHEMBL201': [ // EGFR
      {
        molecule_chembl_id: 'CHEMBL123',
        canonical_smiles: 'CCO',
        pchembl_value: 7.2,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL456',
        canonical_smiles: 'CC(C)O',
        pchembl_value: 6.8,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL789',
        canonical_smiles: 'C1CCCCC1',
        pchembl_value: 6.5,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      }
    ],
    'CHEMBL203': [ // ERBB2
      {
        molecule_chembl_id: 'CHEMBL111',
        canonical_smiles: 'CCN(CC)CC',
        pchembl_value: 7.0,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL222',
        canonical_smiles: 'CC(C)(C)O',
        pchembl_value: 6.9,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL333',
        canonical_smiles: 'C1=CC=CC=C1',
        pchembl_value: 6.7,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      }
    ],
    'CHEMBL197': [ // FGFR1
      {
        molecule_chembl_id: 'CHEMBL444',
        canonical_smiles: 'CC(=O)O',
        pchembl_value: 6.8,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL555',
        canonical_smiles: 'CC(C)(C)C',
        pchembl_value: 6.6,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      },
      {
        molecule_chembl_id: 'CHEMBL666',
        canonical_smiles: 'C1=CC=C(C=C1)O',
        pchembl_value: 6.4,
        assay_type: 'IC50',
        target_organism: 'Homo sapiens'
      }
    ],
    'CHEMBL2363': [ // InhA (TB)
      {
        molecule_chembl_id: 'CHEMBL777',
        canonical_smiles: 'CC(C)CO',
        pchembl_value: 7.1,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL888',
        canonical_smiles: 'CC(C)(C)CO',
        pchembl_value: 6.9,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL999',
        canonical_smiles: 'C1=CC=C(C=C1)CO',
        pchembl_value: 6.7,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      }
    ],
    'CHEMBL2364': [ // KatG (TB)
      {
        molecule_chembl_id: 'CHEMBL101',
        canonical_smiles: 'CC(C)(C)C(=O)O',
        pchembl_value: 6.8,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL202',
        canonical_smiles: 'CC(C)C(=O)O',
        pchembl_value: 6.6,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL303',
        canonical_smiles: 'C1=CC=C(C=C1)C(=O)O',
        pchembl_value: 6.4,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      }
    ],
    'CHEMBL2365': [ // RpoB (TB)
      {
        molecule_chembl_id: 'CHEMBL404',
        canonical_smiles: 'CC(C)(C)C(=O)N',
        pchembl_value: 7.0,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL505',
        canonical_smiles: 'CC(C)C(=O)N',
        pchembl_value: 6.8,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      },
      {
        molecule_chembl_id: 'CHEMBL606',
        canonical_smiles: 'C1=CC=C(C=C1)C(=O)N',
        pchembl_value: 6.6,
        assay_type: 'IC50',
        target_organism: 'Mycobacterium tuberculosis'
      }
    ]
  };
  
  const activities = mockActivities[target] || [];
  
  res.json({
    target: target,
    compound: compound,
    activities: activities,
    count: activities.length
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`ChEMBL MCP Server running on port ${port}`);
});
