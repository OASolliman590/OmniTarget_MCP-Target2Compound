const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock Reactome MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'Reactome MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/species', (req, res) => {
  res.json({
    species: [
      { id: '9606', name: 'Homo sapiens' },
      { id: '10090', name: 'Mus musculus' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: 'R-HSA-123456', name: 'Mock Pathway', species: 'Homo sapiens' }
    ]
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`Reactome MCP Server running on port ${port}`);
});
