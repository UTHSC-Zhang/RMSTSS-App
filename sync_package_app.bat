@echo off
setlocal enabledelayedexpansion

:: Path to Git Bash
set "BASH=C:\Program Files\Git\bin\bash.exe"

:: Relative paths
set "SRC=."
set "DEST=../UTHSC-package/inst/shiny_app"
set "REPO=../UTHSC-package"
set "BRANCH=main"

echo ============================================
echo   Starting Sync and Git Push from UTHSC-App
echo ============================================

:: Fetch and pull latest from remote
echo.
echo 🔄 Pulling latest changes...
"%BASH%" -c "cd '%REPO%' && git fetch origin && git pull origin %BRANCH%"

:: Sync .R and .Rmd files while preserving structure
echo.
echo 🔁 Syncing .R and .Rmd files to shiny_app...
"%BASH%" -c "rsync -av --include '*/' --include '*.R' --include '*.Rmd' --exclude '*' '%SRC%/' '%DEST%/'"

:: Add and commit
echo.
echo ✅ Staging and committing changes...
"%BASH%" -c "cd '%REPO%' && git add inst/shiny_app && git commit -m 'Sync .R and .Rmd files from UTHSC-App on $(date +\"%%Y-%%m-%%d %%H:%%M:%%S\")' || echo 'No changes to commit.'"

:: Push to GitHub
echo.
echo 🚀 Pushing to GitHub...
"%BASH%" -c "cd '%REPO%' && git push origin %BRANCH%"

:: Done
echo.
echo ✅ Done. All operations complete.

pause
