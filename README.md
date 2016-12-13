# ObservationTools
A set of tools to plan astronomical observations

People are more than welcome to do pull requests, open issues, give suggestions, etc.
If you do not have a github user (or don't want to use github for some obscure reason), I can be contacted here: `daniel.andreasen@astro.up.pt`.

# Installation

## With git
If you have git installed (or if you want to install git and use it for the first time), then the tools can be installed with the following few commands in the terminal

    git clone https://github.com/iastro-pt/ObservationTools
    cd ObservationTools
    pip install -r requirements.txt  # You may need to use sudo here


## Without git
If you do not have git installed, you can just download the entire directory [here](https://github.com/iastro-pt/ObservationTools/archive/master.zip).

    unzip ObservationTools-master.zip
    cd ObservationTools-master
    pip install -r requirements.txt  # You may need to use sudo here

## Optional
There is the ability for rv.py to use the c-compiled [ajplant]{https://github.com/andres-jordan/ajplanet.git} module for improved speed.
If ajplanet cannot be installed don't worry because it can also run without ajplanet.

    git clone https://github.com/andres-jordan/ajplanet.git
    cd ajplnaet
    make
    cd Python
    python setup.py build
    python setup.py install

# Updates
If you want to update your tools and installed it with `git`, simply change the directory to this folder and do a `git pull`.
If you don't used git, you have to do the installation again as described above.
