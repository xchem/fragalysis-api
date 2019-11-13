__doc__ = """ 
This is for test purposes only. 

After running this file, query localhost:5000/ 
In a unrelated notebook test it as follows… where directory variable point to the folder you want to test. 

asdasdasdsd

import requests 
import os 

directory = 'data/ATAD2/' 
files = [('files', open(f'{directory}/{file}', 'rb')) for file in os.listdir(directory)] 
r = requests.post('http://localhost:5000/', files=files) 
print(r.json()) 

What just happened is that you sent web request on port 500 where flask is running sending the files. 
"""


from flask import Flask, request
from align import Align
#from .code.main import bla bla
import os

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = '.'

def is_valid_filename(filename):
    """

    :param filename: the user uploaded file
    :type filename: str
    :return: Boolean
    """
    return True

def safe_filename():
    """
    Unique and pymol safe?
    :return: str of file.
    """
    return 'test.pdb'



@app.route('/', methods=['GET','POST'])
def main_route():
    if request.method == 'POST':
        uploaded_files = request.files.getlist("files")
        a_dir = 'test/'

        if not os.path.exists(a_dir):
            os.makedirs(a_dir)

        for file in uploaded_files:
            if not is_valid_filename(file.filename):
                return {'status': f'{file.filename} is not a valid file'}

            filename = safe_filename()
            print(filename)
            print(file.filename)
            file.save(os.path.join(a_dir, filename))

        struc = Align(a_dir, pdb_ref='')
        struc.align()
        return {'status': 'success'}
    else:
        return __doc__

if __name__ == '__main__':
    app.run()
