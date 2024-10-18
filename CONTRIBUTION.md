# Contribution Guidelines

## Setting up user information

Please, have **Git** set up with consistent user information before commiting. Preferably, provide your real name and a working email address.

**Example**:

```bash
git config --global user.name "Your Name Comes Here"
git config --global user.email you@yourdomain.example.com
```

## Make sure that your branch contains clean commits

- Follow the common sense guidelines for writing good commit messages (see below).
- Make separate commits for separate changes. If you cannot describe what the commit does in one sentence, it is probably a mix of changes and should be separated into several commits.
- Do not merge `master` into your branch. Instead use `git rebase`. If you need to resolve merge conflicts or include the latest changes.

## Check your coding style

- Make sure your contributions and changes follow the coding and indentation style of the code surrounding your changes.
- Do not commit commented-out code or files that are no longer needed. Remove the code or the files unless there is a good reason to keep it.

## Guidelines for good commit messages

1. Separate subject from body with a blank line
2. Use the imperative mood in the subject line ("Fix", "Add", "Change" instead of "Fixed", "Added", "Changed")
3. Limit the subject line to 50 characters
4. Reference an issue at the end of a subject line
5. Do not end the subject line with a period
6. Wrap the body at 72 characters
7. Use the body to explain what and why vs. how

Bad:

```bash
some changes fixing this and that...
```

Good:

```bash
fix broken link in reports #123

As foo changed their internal data structure in the last release bar, we need to update our external links accordingly.
```

## Helpful Sources

- https://www.atlassian.com/git/tutorials

## Local testing of Bakta

For local tests and debugging of your contributions, you can create a dev Conda/Mamba environment pointing to local source code files:

1. Clone the repository: `git clone https://github.com/oschwengers/bakta.git`
2. Enter the repo: `cd bakta`
3. Set the conda environment with the dependencies: `mamba create -n bakta-dev -c conda-forge -c bioconda bakta`
4. Activate the environment: `mamba activate bakta-dev`
5. Install via pip pointing to the local source code: `python -m pip install -e .`
