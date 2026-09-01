#!/bin/bash
# Publish web/ to the gh-pages branch as a SINGLE commit.
#
# The bundle is ~17 MB of binaries that are rewritten wholesale on every
# rebuild.  Committing it normally would grow history by that much each time
# and never shrink, so the branch is rebuilt from scratch and force-pushed:
# it always has exactly one commit, and the old blobs fall out of reach.
#
#   ./deploy.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
BRANCH=gh-pages
cd "$HERE"
git rev-parse --git-dir > /dev/null 2>&1 || { echo "not a git repository"; exit 1; }
git remote get-url origin > /dev/null 2>&1 || { echo "no 'origin' remote"; exit 1; }
[ -f web/index.html ] || { echo "no web/index.html"; exit 1; }
[ -f web/.nojekyll ] || { echo "web/.nojekyll missing (Pages would drop _flux etc.)"; exit 1; }

# Build on a throwaway branch and push it to gh-pages by refspec.  Checking
# out `gh-pages` itself would collide with the local branch left by the last
# deploy, and deleting that branch each time is one more thing to get wrong.
TMP=$(mktemp -d)
SCRATCH="deploy-$$"
git worktree add --detach -q "$TMP"
(
  cd "$TMP"
  git checkout -q --orphan "$SCRATCH"
  git rm -rqf . 2>/dev/null || true
  # -R and the trailing dot copy the CONTENTS of web/, so the site sits at the
  # root of the branch and .nojekyll lands where Pages looks for it
  cp -R "$HERE"/web/. .
  git add -A
  git commit -q -m "site $(date -u '+%Y-%m-%d %H:%M UTC')"
  git push -f -q origin "HEAD:$BRANCH"
)
git worktree remove --force "$TMP"
git branch -qD "$SCRATCH" 2>/dev/null || true
echo "pushed $BRANCH  ($(du -sh "$HERE"/web | cut -f1))"
