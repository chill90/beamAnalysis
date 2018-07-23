#Run this script from the top-most directory of the tree you want to clean.
import os

targets = ('~', '.pyc', '.ds_store')

for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if any(target in file.lower() for target in targets):
            print "Deleting file: " + file
            os.remove(os.path.join(root, file))
