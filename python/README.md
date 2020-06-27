# Usage

# Dependencies
Dependencies are listed pyproject.toml and poetry.lock contains the full list of
dependencies and sub-dependencies. To easily add dependencies and resolve dependency conflicts,
we use poetry https://python-poetry.org/docs/.

# Developer notes
To export poetry.lock content to requirements files, use:

´poetry export -f requirements.txt > requirements.txt --without-hashes´ 