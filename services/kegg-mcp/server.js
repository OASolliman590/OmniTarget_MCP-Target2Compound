const express = require('express');
const app = express();
const port = 3000;

app.use(express.json());

// Mock KEGG MCP endpoints
app.get('/info', (req, res) => {
  res.json({
    name: 'KEGG MCP Server',
    version: '1.0.0',
    status: 'running'
  });
});

app.get('/search/pathway', (req, res) => {
  res.json({
    query: req.query.q || 'depression',
    results: [
      { id: 'hsa05032', name: 'Morphine addiction', species: 'hsa' },
      { id: 'hsa05030', name: 'Cocaine addiction', species: 'hsa' }
    ]
  });
});

app.get('/pathway/:id', (req, res) => {
  res.json({
    id: req.params.id,
    name: 'Mock Pathway',
    genes: ['GENE1', 'GENE2', 'GENE3']
  });
});

app.listen(port, '0.0.0.0', () => {
  console.log(`KEGG MCP Server running on port ${port}`);
});
