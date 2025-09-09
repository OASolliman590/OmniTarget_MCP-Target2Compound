const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock ProteinAtlas MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'ProteinAtlas MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/genes', (req, res) => {
  res.json({
    genes: [
      { id: 'ENSG00000000003', name: 'TSPAN6', description: 'Tetraspanin 6' },
      { id: 'ENSG00000000005', name: 'TNMD', description: 'Tenomodulin' }
    ]
  });
});

app.get('/search/query', (req, res) => {
  res.json({
    query: req.query.q || 'test',
    results: [
      { id: 'ENSG00000000003', name: 'TSPAN6', tissue: 'Brain' }
    ]
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`ProteinAtlas MCP Server running on port ${port}`);
});
