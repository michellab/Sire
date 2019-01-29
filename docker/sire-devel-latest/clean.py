import os

py_saved_space = 0
print("Removing cached .pyc and .pyo files...")
for root, dirs, files in os.walk("/home/sireuser/sire.app"):
    for file in files:
        if file.endswith(".pyc") or file.endswith(".pyo"):
            filename = os.path.join(root,file)
            py_saved_space += os.path.getsize(filename)
            os.remove(filename)

print("...files removed. Saved %d MB of space" % (py_saved_space/(1024*1024)))

