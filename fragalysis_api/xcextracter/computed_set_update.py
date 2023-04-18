import os
import sys
import requests
import time
import threading
import _thread as thread

REQ_URL = 'https://fragalysis.diamond.ac.uk/viewer/upload_cset/'

# The keycloak credentials (sensitive).
# These are required if you're going to use the authenticated API.
# They will be asserted when the appropriate functions are used.
KEYCLOAK_USERNAME = os.environ.get('KEYCLOAK_USERNAME')
KEYCLOAK_PASSWORD = os.environ.get('KEYCLOAK_PASSWORD')
KEYCLOAK_CLIENT_SECRET = os.environ.get('KEYCLOAK_CLIENT_SECRET')


def get_csrf(REQ_URL):
    """Get a csrf token from the request url to authenticate further requests
    Parameters
    ----------
    REQ_URL string
        The URL that you want to make a request against after getting the token
    
    Returns
    -------
    csrftoken
        csrf token to use for further requests against the same URL
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


def get_keycloak_access_token(*,
                              keycloak_url,
                              keycloak_realm,
                              keycloak_client_id):
    """Gets an 'access token' from Keycloak.
    If successful we'll get the token (a big long string).
    
    Each access tokens should last for a number of minutes depending on the
    keycloak configuration.
    
    To use this function you must have set the corresponding
    environment variables KEYCLOAK_USERNAME, KEYCLOAK_PASSWORD,
    and KEYCLOAK_CLIENT_SECRET
    
    Parameters
    ----------
    keycloak_url string
        The keycloak authentication URL (ending in '/auth')
    keycloak_realm string
        The keycloak realm the stack belongs to
    keycloak_client_id string
        The client ID the stack is known by in keycloak

    Returns
    -------
    access_token string
        An API access token
    """
    assert KEYCLOAK_USERNAME
    assert KEYCLOAK_PASSWORD
    assert KEYCLOAK_CLIENT_SECRET
    
    print(KEYCLOAK_USERNAME)
    print(KEYCLOAK_PASSWORD)
    print(KEYCLOAK_CLIENT_SECRET)
    
    realm_url = f'{keycloak_url}/realms/{keycloak_realm}'
    url = f'{realm_url}/protocol/openid-connect/token'
    data = (
        f'client_id={keycloak_client_id}'
        f'&grant_type=password'
        f'&username={KEYCLOAK_USERNAME}'
        f'&password={KEYCLOAK_PASSWORD}'
        f'&client_secret={KEYCLOAK_CLIENT_SECRET}'
    )
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    resp = requests.post(url, headers=headers, data=data, timeout=4.0)
    assert resp.status_code == 200
    assert 'access_token' in resp.json()
    return resp.json()['access_token']

def upload_cset(REQ_URL, access_token, target_name, sdf_path, upload_key=None, pdb_zip_path=None):
    """Upload a computed set to Fragalysis

    Parameters
    ----------
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    access_token: str
        a valid OIDC/Keycloak access token, inserted into the request header as a bearer token
    target_name: str
        the name of the target in Fragalysis that the computed set is for
    sdf_path: str
        path to the sdf file to upload
    upload_key: str
        upload key, not currently turned on, so can be any value, but not blank or null (optional)
    pdb_zip_path: str
        path to the zip file of pdb's to upload (optional)
    """
    assert access_token

    print(f'Uploading new computed set files....')
    
    csrf_token = get_csrf(REQ_URL)
    
    payload = {'target_name': target_name,
                'submit_choice': "U",
                'upload_key': upload_key,
                'update_set': None}

    files = [
        ('sdf_file', open(sdf_path,'rb')),
        
    ]
    
    if pdb_zip_path:
        files.append(('pdb_zip', open(pdb_zip_path,'rb')))

    headers = {'X-CSRFToken': csrf_token,
              'Cookie': f'csrftoken={csrf_token}',
              'Authorization': f'Bearer {access_token}'}
    
    response = requests.request("POST", REQ_URL, headers=headers, data=payload, files=files)
    
    taskurl = get_taskidurl(response=response, REQ_URL=REQ_URL)
    
    if not taskurl:
        raise Exception(f'Something went wrong with the upload request! \
                Please try again or email warren.thompson@diamond.ac.uk for help.\
                Response: {response.text}') 
    
    return taskurl


def update_cset(REQ_URL, access_token, target_name, sdf_path, cset_name, upload_key=None, pdb_zip_path=None, add=False):
    """Send data to <root_url>/viewer/upload_cset/ to overwrite an existing computed set, or to 
    <root_url>/viewer/update_cset/ to add new molecules without deleting the old ones.

    Parameters
    ----------
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    access_token: str
        a valid OIDC/Keycloak access token, inserted into the request header as a bearer token
    target_name: str
        the name of the target in Fragalysis that the computed set is for
    cset_name: str
        the name of the computed set you want to update 
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file). Leave blank if you are adding
        a set for the first time
    sdf_path: str
        path to the sdf file to upload
    upload_key: str
        upload key, not currently turned on, so can be any value, but not blank or null (optional)
    pdb_zip_path: str
        path to the zip file of pdb's to upload (optional)
    add: bool
        set to True if updating a computed set without overwriting it completely (for <root_url>/viewer/update_cset/)

    Returns
    -------
    taskurl: str
        The URL of the Celery task
    """
    assert access_token

    print(f'Submitting files to update {cset_name}...')
    
    csrf_token = get_csrf(REQ_URL)
        
    if not add:
        payload = {'target_name': target_name,
                   'submit_choice': "U",
                   'upload_key': upload_key,
                   'update_set': cset_name}
    else:
        payload = {'target_name': target_name,
                   'update_set': cset_name}

    files = [
        ('sdf_file', open(sdf_path,'rb')),
        
    ]
    
    if pdb_zip_path:
        files.append(('pdb_zip', open(pdb_zip_path,'rb')))

    headers = {'X-CSRFToken': csrf_token,
              'Cookie': f'csrftoken={csrf_token}',
              'Authorization': f'Bearer {access_token}'}
    
    response = requests.request("POST", REQ_URL, headers=headers, data=payload, files=files)
    
    taskurl = get_taskidurl(response=response, REQ_URL=REQ_URL)
    
    if not taskurl:
        raise Exception(f'Something went wrong with the update request! \
                Please try again or email warren.thompson@diamond.ac.uk for help.\
                Response: {response.text}') 
    
    return taskurl


def delete_cset(REQ_URL, access_token, target_name, cset_name, upload_key=None):
    """Delete uploaded computer set

    Parameters
    ----------
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    access_token: str
        a valid OIDC/Keycloak access token, inserted into the request header as a bearer token
    target_name: str
        the name of the target in Fragalysis that the computed set is for
    cset_name: str
        the name of the computed set you want to delete 
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file). Leave blank if you are adding
        a set for the first time
    upload_key: str
        upload key, not currently turned on, so can be any value, but not blank or null (optional)
    """
    assert access_token
    
    print(f'Deleting {cset_name}...')
    
    csrf_token = get_csrf(REQ_URL)
        
    payload = {'target_name': target_name,
                'submit_choice': "D",
                'upload_key': upload_key,
                'update_set': cset_name}

    headers = {'X-CSRFToken': csrf_token,
              'Cookie': f'csrftoken={csrf_token}',
              'Authorization': f'Bearer {access_token}'}
    
    response = requests.request("POST", REQ_URL, headers=headers, data=payload)
    
    check_url_csetdeleted(cset_name=cset_name, response=response)

def check_url_csetdeleted(cset_name: str, response: str):
    """Checks if response yields delete success

    Parameters
    ----------
    cset_name: str
        the name of the computed set you want to delete 
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file). Leave blank if you are adding
        a set for the first time
    response: str
        The response from the update cset call
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    """
    lines = response.text.split('\n')
    delete_check = None
    for l in lines:
        if f'<p style="color:blue;">Compound set &quot;{cset_name}&quot; deleted</p>' in l:
            delete_check = True
            print(f'Computed set: {cset_name} deleted')
            break
    if not delete_check:
        raise Exception(f'Something went wrong with the delete request! \
                        Please try again or email warren.thompson@diamond.ac.uk for help.\
                        Response: {response.text}')             

def get_taskidurl(response: str, REQ_URL: str):
    """Checks if response yields a Celery task ID and returns the URL

    Parameters
    ----------
    response: str
        The response from the update cset call
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    """
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
    return taskurl


def quit_function(fn_name):
    """Quit a function and return an error
    Parameters
    ----------
    fn_name:
        name of function to apply to
    """
    # print to stderr, unbuffered in Python 2.
    print('{0} took too long. The task has probably not worked, but is left in a PENDING state. \
    Please try again or email rachael.skyner@diamond.ac.uk for help.'.format(fn_name), file=sys.stderr)
    sys.stderr.flush() # Python 3 stderr is likely buffered.
    thread.interrupt_main() # raises KeyboardInterrupt
    

def exit_after(s):
    '''
    use as decorator to exit process if function takes longer than s seconds
    
    Parameters
    ----------
    s: int
        number of seconds to exit after
    '''
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
    """Check a task url to get it's status. Will return SUCCESS or FAILED, or timeout after 10 minutes (600s)
       if the task is still pending
    
    Parameters
    ----------
    taskurl: str
        URL to ping
        
    Returns
    -------
    status: str
        SUCCESS or FAILURE
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


# EXAMPLES:
# =========
# ---- NB: The major difference here is the REQ_URL. For new data, or to overwrite data, use https://fragalysis.diamond.ac.uk/viewer/upload_cset/.
#          To add new molecules to an existing set, use https://fragalysis.diamond.ac.uk/viewer/update_cset/. ----
#
# You must set the following environment variables
# to avoid assertions when you run the code:
#
# - KEYCLOAK_USERNAME (a valid user for the chosen stack)
# - KEYCLOAK_PASSWORD (the user's password)
# - KEYCLOAK_CLIENT_SECRET (this is a stack-specific uuid4-like value)
#
# i.e.:
#
#   export KEYCLOAK_USERNAME=someone
#   export KEYCLOAK_PASSWORD=someone1234
#   export KEYCLOAK_CLIENT_SECRET=00000000-0000-0000-0000-000000000000
#
# to get an OIDC/Keycloak access token:
# -------------------------------------
# Here we get a token for the staging stack...
#
# access_token = get_keycloak_access_token(
#     keycloak_url = 'https://keycloak.xchem.diamond.ac.uk/auth',
#     keycloak_realm = 'xchem',
#     keycloak_client_id = 'fragalysis-xchem')
# print(f'access_token="{access_token}"')
#
# to upload a new cset:
# ------------------------------
# taskurl = upload_cset(REQ_URL='https://fragalysis.xchem.diamond.ac.uk/viewer/upload_cset/',
#                       access_token=access_token,
#                       target_name='Mpro',
#                       upload_key='1',
#                       sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28.sdf',
#                       pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs.zip')
# task_response, json_results = get_task_response(taskurl)
# print(task_response)
# print(json_results)
#
# to overwrite an existing cset:
# ------------------------------
# taskurl = update_cset(REQ_URL='https://fragalysis.xchem.diamond.ac.uk/viewer/upload_cset/',
#                       access_token=access_token,
#                       target_name='Mpro',
#                       upload_key='1',
#                       cset_name='WT-xCOS3-ThreeHop',
#                       sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28.sdf',
#                       pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs.zip')
# task_response, json_results = get_task_response(taskurl)
# print(task_response)
# print(json_results)
#
# to update an existing cset:
# ---------------------------
# taskurl = update_cset(REQ_URL='https://fragalysis.xchem.diamond.ac.uk/viewer/update_cset/',
#                       access_token=access_token,
#                       target_name='Mpro',
#                       cset_name='WT-xCOS3-ThreeHop',
#                       sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28 copy.sdf',
#                       pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs copy.zip',
#                       add=True)
# task_response, json_results = get_task_response(taskurl)
# print(task_response)
# print(json_results)
#
# to delete an existing cset:
# ---------------------------
# delete_cset(REQ_URL='https://fragalysis.xchem.diamond.ac.uk/viewer/update_cset/',
#             access_token=access_token,
#             target_name='Mpro',
#             cset_name='WT-xCOS3-ThreeHop',
#             upload_key='1')
# task_response, json_results = get_task_response(taskurl)
# print(task_response)
# print(json_results)