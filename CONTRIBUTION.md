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
4. Reference an issue in the beginning of a subject line
5. Do not end the subject line with a period
6. Wrap the body at 72 characters
7. Use the body to explain what and why vs. how

Bad:
<pre>
Some changes fixing this and that...
</pre>

Good:
<pre>
#123 Fix broken link in ABR reports

As CARD changed some of their website structure we need to update link templates for antibiotic resistance (ABR) models.
</pre>

# Helpful Sources
- https://www.atlassian.com/git/tutorials

