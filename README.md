client side:
pip install git+https://github.com/mattpre/zsoil_tools.git

for development:
got to directory where source code should be downloaded
git clone https://github.com/mattpre/zsoil_tools.git
cd zsoil_tools
python setup.py develop (might have to provide full path to python.exe)

aide-memoire MP:
python setup.py bdist_egg
git add .
git commit -m"some message"
git push

