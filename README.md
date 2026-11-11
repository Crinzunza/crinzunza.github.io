# crinzunza.github.io

Personal academic website of Cristian Inzunza, built with plain Jekyll (no theme) and served by GitHub Pages. English pages live at the root, Spanish pages under `/es/`.

## Where to edit

| To change | Edit |
|---|---|
| Publications (both languages, home page, research and Fondecyt pages) | `_data/publications.yml` |
| News | `_data/news.yml` |
| Talks | `_data/talks.yml` |
| Teaching | `_data/teaching.yml` |
| Research lines (texts, captions, figures) | `_data/research.yml` |
| Biography | `index.html` (EN) and `es/index.html` (ES) |
| Fondecyt project text | `fondecyt.html` and `es/fondecyt.html` |
| Name, email, affiliation, profile links | `_config.yml` |
| Menu labels and interface strings | `_data/i18n.yml` |
| CV | replace `files/cv.pdf` |
| Colours and typography | `assets/css/site.css` (tokens at the top) |

### Adding a publication

Add an entry at the top of `_data/publications.yml`:

```yaml
- title: "..."
  authors: "A. Author, C. Inzunza"
  venue: "Journal name, volume, pages"
  year: 2027
  status: published        # or: review
  topics: [poro, apost]    # ids from _data/research.yml
  doi: "10.xxxx/xxxxx"     # optional; also arxiv, preprint, code
  fondecyt: true           # optional; lists it on the Fondecyt page
```

When a preprint is accepted, change `status` to `published` and add the `doi`.

### Using a real simulation figure

Put the image in `assets/img/` and add `image: /assets/img/name.png` to the research line in `_data/research.yml`. It replaces the generated mesh drawing.

### Figures

The mesh drawings (`_includes/hero-mesh.svg`, `_includes/thumb-*.svg`) were generated with a Delaunay triangulation script (NumPy + SciPy). Their colours come from CSS variables, so they follow light and dark mode.

## Local preview

```bash
bundle install
bundle exec jekyll serve
```

Then open http://localhost:4000.
