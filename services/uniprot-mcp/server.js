const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock UniProt MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'UniProt MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/status', (req, res) => {
  res.json({
    status: 'running',
    name: 'UniProt MCP Server',
    description: 'Mock UniProt database server'
  });
});

app.get('/genes', (req, res) => {
  res.json({
    genes: [
      { id: 'P12345', name: 'BRCA1', description: 'Breast cancer type 1 susceptibility protein' },
      { id: 'P38398', name: 'BRCA2', description: 'Breast cancer type 2 susceptibility protein' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: 'P12345', name: 'BRCA1', organism: 'Homo sapiens' }
    ]
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`UniProt MCP Server running on port ${port}`);
});
