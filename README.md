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

## Initial setup:
- install Quarto:

    ```bash
    brew install quarto
    ```

- install `blogdown`

    ```bash
    Rscript -e "install.packages('blogdown', repos='https://cran.rstudio.com/')
    ```

- install `jupyter`:

    ```bash
    pip3 install jupyter jupyter-cache
    ```

## Runtime dependencies

- For R, I usually just use the system R
- For Python,
    Set the appropriate kernel by registering it with Jupyter

    ```python
    python -m ipykernel install --name <<KERNEL_NAME>> --display-name <<DISPLAY_NAME>> --user
    ```
    
    Add the kernel to the yaml frontmatter.
    
    ```yaml
    jupyter: <<KERNEL_NAME>>
    ```
    
## Converting notebooks to markdown

1. Source scripts/build.py to knit .Rmd and .qmd documents in _source into .md documents in _posts
2. `bundle exec jekyll serve` to create static site in the _site directory
3. Fix LaTeX errors (using packages which aren't recognized by MathJaX may throw errors in build.R, while rendering mistakes are also common)

# Site Template

- Forked from [Minimal Mistakes](https://github.com/mmistakes/minimal-mistakes):
- Web site design derived from [Lanyon](https://github.com/poole/lanyon) by Mark Otto.
- Open sourced under the [MIT license](LICENSE.md).
