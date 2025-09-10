const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock PDB MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'PDB MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/status', (req, res) => {
  res.json({
    status: 'running',
    name: 'PDB MCP Server',
    description: 'Mock PDB database server'
  });
});

app.get('/genes', (req, res) => {
  res.json({
    genes: [
      { id: '1A2B', name: 'Mock Protein', description: 'Mock protein structure' },
      { id: '2C3D', name: 'Mock Enzyme', description: 'Mock enzyme structure' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: '1A2B', name: 'Mock Protein', resolution: 2.5 }
    ]
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`PDB MCP Server running on port ${port}`);
});
