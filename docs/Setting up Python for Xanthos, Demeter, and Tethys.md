# Setting up Python for Xanthos, Tethys, and Demeter

Last revised: 10/5/2020

## Installations

  1. Follow all the installations in installations doc, these are just the relevant ones for Python.
  2. Install PyCharm - Professional version: https://www.jetbrains.com/pycharm/. Get license from https://www.jetbrains.com/community/education/\#students and use it to activate.
  3. Install Python 64-bit (Python 3.8 is needed for Xanthos and Tethys, Python 2.7 is needed for Demeter).
     * Save to Users\yourusername\AppData\Local\Programs\Python. Here we will have 2 folders "Python38 and Python27"  
     * Note: AppData is a hidden folder so you may need to type it out (rather than clicking in the window)  
     * Make sure the python.exe is saved  
  4. Install Git LFS a. In Git Bash, type ```git lfs install``` to initialize Git LFS.
  
## Clone Xanthos, Tethys, and Demeter

  5. Clone Xanthos, Tethys and Demeter. Recommend to clone these into a folder (e.g. called "INFEWS"" for the rest of this document) hard drive with a lot of space.

## Setting up new Python environment in PyCharm

  6. Open this folder ("INFEWS"") in PyCharm as a project.
  7. Go to File \> Settings \> Project: INFEWS \> Python Interpreter.
  8. Click the Settings icon \> Add.
  9. In Location: enter D:\INFEWS\venv38 or \venv27, depending on if you are installing Python 3.8 or Python 2.7 respectively. We want to have 2 separate environments for these and be able to switch between them easily.
  10. In Base interpreter: navigate to the respective python.exe in the Python38 or Python27 folder.
  11. Click "OK".
  12. Under packages, make sure "pip"" and "setuptools"" are installed.

## Adding Xanthos, Tethys and Demeter to Interpreter Paths

  13. We also want to add the xanthos, tethys, and demeter folder path to Interpreter Paths because part of corresponding packages are installed within those folders. Go to File \> Settings \> Project: INFEWS.
  14. Click the Settings icon \> Show All.
  15. Click the version of Python we are adding to. For tethys and xanthos, it is Python 3.8 and for demeter it is Python 2.7.
  16. Click the lowest button on the right side of the window (Show paths for the selected interpreter).
  17. Click the + icon.
  18. Navigate to the tethys, xanthos or demeter folder within the INFEWS folder and Click "OK".

## Running the setup.py install

  19. Go to Terminal and navigate to one of the cloned folders (tethys, xanthos, etc).
  20. Run ```python setup.py install```.
  21. It should add Build and other folders to your cloned folder.

**Notes for trouble shooting**
  * For setting up Demeter, may run into issues because not all the packages are installed.
  * If there is an error within a joblib subfolder ("Missing package: setuptools not found. Demeter requires this to install. Please install setuptools and retry."). It's actually the case that other packages are not installed, even though setuptools is there. This is due to using Python 2.7 which has errors installing and updating packages. Need to manually install packages.

  22. Install the following packages.
  
      ```pip install joblib```
      
  23. Using the --use-feature=2020-resolver, install remaining packages needed (e.g. configobj, matplotlib, pyshp, scipy).
  
      ```pip install matplotlib --use-feature=2020-resolver```

## Change file paths

  24. For Xanthos, we need to change the file paths in the "pm\_abcd\_mrtm.ini" file.
      * Add the RootDir (give the full path of your xanthos/example folder, e.g. D:/INFEWS/xanthos/emaple).
      * Change other file paths accordingly to point to the correct path for xanthos/example (or whatever script you are running).
  25. Similarly for Demeter, we need to update file paths in the "config.ini" file.

***Try running the example.py script!***
