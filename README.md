# NC Interpolation — Interpolation Calculator

A web-based interpolation calculator with step-by-step solutions, interactive plots, and multiple method support.

## Methods

- **Lagrange Interpolation** — For unequally spaced data
- **Newton's Divided Difference** — For unequally spaced data
- **Newton's Forward Difference** — For equally spaced data
- **Newton's Backward Difference** — For equally spaced data
- **Stirling's Formula** — Central difference interpolation

## Features

- Input data points (x, y pairs) and value to interpolate
- Auto-detect equally/unequally spaced data
- Step-by-step solutions with difference tables
- Interactive polynomial curve plotting on canvas
- Compare mode for multiple methods side-by-side
- Dark/light theme with system preference detection
- Calculation history with localStorage persistence
- Share links for reproducible calculations
- LaTeX export for academic use
- Print-friendly output
- Responsive design for mobile and desktop
- Full SEO with meta tags, OG, Twitter cards, JSON-LD

## Tech Stack

- Vanilla HTML, CSS, JavaScript
- Math.js from CDN (for future extensions)
- Docker + nginx for deployment

## Development

Open `index.html` in a browser. No build step required.

## Deployment

```bash
docker compose up -d --build
```

## License

MIT
