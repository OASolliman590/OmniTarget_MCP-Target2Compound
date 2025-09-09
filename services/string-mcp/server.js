const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock STRING MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'STRING MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/version', (req, res) => {
  res.json({
    version: '1.0.0',
    name: 'STRING MCP Server',
    description: 'Mock STRING database server'
  });
});

app.get('/genes', (req, res) => {
  res.json({
    genes: [
      { id: '9606.ENSP00000000233', name: 'ARF5', description: 'ADP ribosylation factor 5' },
      { id: '9606.ENSP00000000412', name: 'M6PR', description: 'Mannose-6-phosphate receptor' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: '9606.ENSP00000000233', name: 'ARF5', score: 0.95 }
    ]
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`STRING MCP Server running on port ${port}`);
});
