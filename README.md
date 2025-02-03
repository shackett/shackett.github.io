# shackett.org

# One Time Setup

```bash
# install jekyll
gem install jekyll bundler

# install dependencies from Gemfile.lock
bundle install

# serve site locally
bundle exec jekyll serve
```

# Adding New Content

1. Source build.R to knit .Rmd documents in _source into .md documents in _posts
2. `bundle exec jekyll serve` to create static site in the _site directory
3. Fix LaTeX errors (using packages which aren't recognized by MathJaX may throw errors in build.R, while rendering mistakes are also common)

# Site Template

- Forked from [Minimal Mistakes](https://github.com/mmistakes/minimal-mistakes):
- Web site design derived from [Lanyon](https://github.com/poole/lanyon) by Mark Otto.
- Open sourced under the [MIT license](LICENSE.md).
