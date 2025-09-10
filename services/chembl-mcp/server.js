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

app.listen(port, '0.0.0.0', () => {
  console.log(`ChEMBL MCP Server running on port ${port}`);
});
