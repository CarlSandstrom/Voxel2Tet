---
name: Commit current changes
description: Create a commit message and commit the current changes in the repository. Use when asked to commit changes.
---

1. Generate a concise commit message summarizing the changes made. 
2. Use `git add .` to stage all changes.
3. Use `git commit -m "<commit message>"` to commit the changes with the generated commit message.
4. Return the commit message used for the commit.
5. If there are no changes to commit, say that explicitly and do not attempt to create a commit.