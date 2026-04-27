# Academic personal site

Static site intended for GitHub Pages at `https://<username>.github.io`.

## Local preview

```sh
python3 -m http.server 8000
# then open http://localhost:8000
```

## Publish

1. Rename this folder to `<your-github-username>.github.io`.
2. Create a public repo on GitHub with the same name.
3. From inside the folder:
   ```sh
   git init
   git add .
   git commit -m "Initial site"
   git branch -M main
   git remote add origin git@github.com:<username>/<username>.github.io.git
   git push -u origin main
   ```
4. The site will be live at `https://<username>.github.io` within a minute or two.

## Files

- `index.html` — page content (about, news, publications, teaching)
- `style.css` — styling
- `assets/` — drop `avatar.jpg`, `cv.pdf`, etc. here
