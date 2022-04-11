if not exist "build\html\NUL" mkdir "build\html"
if not exist "build\doctrees\NUL" mkdir "build\doctrees"
sphinx-build -b html -d build\doctrees -j8 -v source build\html