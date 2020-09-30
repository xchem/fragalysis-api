import sys
import requests
import json
import time
import threading
import _thread as thread

REQ_URL = 'https://fragalysis.diamond.ac.uk/viewer/upload_cset/'

def get_csrf(REQ_URL):
    """! Get a csrf token from the request url to authenticate further requests

    @param REQ_URL The URL that you want to make a request against after getting the token
    @return csrftoken csrf token to use for further requests against the same URL
    """

    client = requests.session()
    # Retrieve the CSRF token first
    client.get(REQ_URL)  # sets cookie
    if 'csrftoken' in client.cookies:
        # Django 1.6 and up
        csrftoken = client.cookies['csrftoken']
    else:
        # older versions
        csrftoken = client.cookies['csrf']
    return csrftoken


def update_cset(REQ_URL, target_name, sdf_path, update_set='None', submit_choice=None, upload_key=None, pdb_zip_path=None, add=False):
    """! Send data to <root_url>/viewer/upload_cset/ to overwrite an existing computed set, or to 
    <root_url>/viewer/update_cset/ to add new molecules without deleting the old ones.

    @param REQ_URL request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    @param target_name the name of the target in Fragalysis that the computed set is for
    @param update_set the name of the computed set you want to update 
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file). Leave blank if you are adding
        a set for the first time
    @param sdf_path str path to the sdf file to upload
    @param submit_choice int 0 for validate, 1 for upload (not required for update - ie. viewer/update_cset)
    @param upload_key upload key, not currently turned on, so can be any value, but not blank or null (optional)
    @param pdb_zip_path path to the zip file of pdb's to upload (optional)
    @param add: bool set to True if updating a computed set without overwriting it completely (for <root_url>/viewer/update_cset/)
    @return taskurl the URL to check for the status of the upload
    """
    print(f'Submitting files to update {update_set}...')
    
    csrf_token = get_csrf(REQ_URL)
    
    if not add:
        payload = {'target_name': target_name,
                   'submit_choice': submit_choice,
                   'upload_key': upload_key,
                   'update_set': update_set}
    else:
        payload = {'target_name': target_name,
                   'update_set': update_set}

    files = [
        ('sdf_file', open(sdf_path,'rb')),
        
    ]
    
    if pdb_zip_path:
        files.append(('pdb_zip', open(pdb_zip_path,'rb')))

    headers = {'X-CSRFToken': csrf_token,
              'Cookie': f'csrftoken={csrf_token}'}
    
    response = requests.request("POST", REQ_URL, headers=headers, data=payload, files=files)
    
    lines = response.text.split('\n')
    taskurl = None
    for l in lines:
        if 'taskUrl = "/viewer/upload_task/' in l:
            taskid = l.split('/')[-2]
            print(f'upload task id: {taskid}')
            taskurl = f'{REQ_URL.replace("/viewer/upload_cset/","/viewer/upload_task/")}{taskid}'
        elif 'taskUrl = "/viewer/update_task/' in l:
            taskid = l.split('/')[-2]
            print(f'upload task id: {taskid}')
            taskurl = f'{REQ_URL.replace("/viewer/update_cset/","/viewer/update_task/")}{taskid}'

            if taskurl:
                break
    if not taskurl:
        raise Exception(f'Something went wrong with the upload/update request! \
                        Please try again or email rachael.skyner@diamond.ac.uk for help.\
                        Response: {response.text}')

    return taskurl


def quit_function(fn_name):
    """! Quit a function and return an error
    
    @param fn_name name of function to apply to
    """
    # print to stderr, unbuffered in Python 2.
    print('{0} took too long. The task has probably not worked, but is left in a PENDING state. \
    Please try again or email rachael.skyner@diamond.ac.uk for help.'.format(fn_name), file=sys.stderr)
    sys.stderr.flush() # Python 3 stderr is likely buffered.
    thread.interrupt_main() # raises KeyboardInterrupt
    

def exit_after(s):
    """! use as decorator to exit process if function takes longer than s seconds
    
    @param s integer number of seconds to exit after
    """
    def outer(fn):
        def inner(*args, **kwargs):
            timer = threading.Timer(s, quit_function, args=[fn.__name__])
            timer.start()
            try:
                result = fn(*args, **kwargs)
            finally:
                timer.cancel()
            return result
        return inner
    return outer


@exit_after(600)
def get_task_response(taskurl):
    """! Check a task url to get it's status. Will return SUCCESS or FAILED,
        or timeout after 10 minutes (600s)if the task is still pending
    
        @param taskurl URL to ping
        
        @return status SUCCESS or FAILURE
        """
        
    print('pinging task to check status...')
    requests.request("GET", taskurl)
    complete=False
    while not complete:
        task_response = requests.request("GET", taskurl)
        if 'upload' in taskurl:
            status = task_response.json()['upload_task_status']
        if 'update' in taskurl:
            status = task_response.json()['update_task_status']
        if status == "SUCCESS":
            complete=True
        if status == "FAILURE":
            complete=True
        time.sleep(5)
    return status, task_response.json()

