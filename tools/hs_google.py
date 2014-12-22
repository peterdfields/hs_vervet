"""
Some tools to interact with the google api
"""
import os, string, itertools
import numpy as np
import pandas as pd
eu = os.path.expanduser
import httplib2
import pprint

from apiclient.discovery import build
from apiclient.http import MediaFileUpload
from apiclient import errors
from oauth2client.client import OAuth2WebServerFlow
from oauth2client.file import Storage

import gspread

def make_credentials():
    # Copy your credentials from the console
    CLIENT_ID = '824765796499-j5a9h3tsrv0f2slnp1961e5bcig15fdj.apps.googleusercontent.com'
    CLIENT_SECRET = 'k9P7ZM92GXvdsEJv3v2TyYkE'

    # Check https://developers.google.com/drive/scopes for all available scopes
    OAUTH_SCOPE = ['https://www.googleapis.com/auth/drive','https://spreadsheets.google.com/feeds', 'https://docs.google.com/feeds']

    # Redirect URI for installed apps
    REDIRECT_URI = 'urn:ietf:wg:oauth:2.0:oob'

    # Path to the file to upload
    FILENAME = 'document.txt'

    # Run through the OAuth flow and retrieve credentials
    flow = OAuth2WebServerFlow(CLIENT_ID, CLIENT_SECRET, OAUTH_SCOPE,
                               redirect_uri=REDIRECT_URI)
    authorize_url = flow.step1_get_authorize_url()
    print 'Go to the following link in your browser: ' + authorize_url
    code = raw_input('Enter verification code: ').strip()
    credentials = flow.step2_exchange(code)
    store = Storage(eu("~/.google/oauth2.dat"))
    store.put(credentials)

def get_credentials():
    OAUTH2_STORAGE = eu('~/.google/oauth2.dat')
    storage = Storage(OAUTH2_STORAGE)
    credentials = storage.get()
    return credentials

def get_drive_service(http,credentials=None):
    # Create an httplib2.Http object and authorize it with our credentials
    if credentials is None:
        credentials = get_credentials()
    drive_service = build('drive', 'v2', http=http)
    return drive_service

def new_document(title, media_body=None, description=None, mimeType='text/plain', parents=None,drive_service=None):
    if drive_service is None:
        drive_service = get_drive_service()
    body = {
      'title': title,
      'description': description,
      'mimeType': mimeType,
       "parents": parents}
    file = drive_service.files().insert(body=body,media_body=media_body).execute()
    return file


def retrieve_all_files(service):
  """Retrieve a list of File resources.

  Args:
    service: Drive API service instance.
  Returns:
    List of File resources.
  """
  result = []
  page_token = None
  while True:
    try:
      param = {}
      if page_token:
        param['pageToken'] = page_token
      files = service.files().list(**param).execute()

      result.extend(files['items'])
      page_token = files.get('nextPageToken')
      if not page_token:
        break
    except errors.HttpError, error:
      print 'An error occurred: %s' % error
      break
  return result

def upload_df(df,spreadsheet,worksheet_title, index=True, header=True):
    """
    upload dataframe  including index and column labels
    to new worksheet within spreadsheet
    """
    if index:
        data = np.hstack([df.index.values.reshape((-1,1)),df.values])
        if header:
            data = np.vstack([np.concatenate([[df.index.name],df.columns.values]),data])
    elif header:
        data = np.vstack([df.columns.values,df.values])
    else:
        data = df.values
    rows = data.shape[0]
    cols = data.shape[1]
    assert cols <= 26, "Not implemented for more than 26 columns."
    #add worksheet
    worksheet = spreadsheet.add_worksheet(title=worksheet_title,rows=str(rows),cols=str(cols))
    cell_list = worksheet.range('A1:{}{}'.format(string.ascii_uppercase[cols-1],rows))
    #add column labels
    for cell, val in itertools.izip(cell_list,data.reshape(-1)):
        cell.value = val
    worksheet.update_cells(cell_list)
    return worksheet


class DriveInteractor(object):
    """
    An object that allows to work with the google drive api
    """
    def __init__(self):
        try:
            self.credentials = get_credentials()
        except:
            raise
        http = httplib2.Http()
        http = self.credentials.authorize(http)
        self.http = http
        self.drive_service = get_drive_service(self.http,credentials=self.credentials)
        self.files = None
        self.gc = None
    
    def refresh_http(self):
        self.credentials.refresh(self.http)
            
    def refresh_drive_service(self):
        self.drive_service = get_drive_service(self.http,credentials=self.credentials)

    def new_document(self, title, media_body=None, description=None, mimeType='text/plain', parents=None):
        return new_document(title, media_body, description, mimeType, parents, drive_service=self.drive_service)


    def get_all_files(self):
        #why not just self.files = service.files().list().execute()
        self.files = retrieve_all_files(self.drive_service)
    

    def id_from_properties(self,**kwa):
        """
        kwas should be file property name and value pairs
        eg.: mimeType='application/vnd.google-apps.folder', title="xpclr"
        """
        if self.files is None:
            self.get_all_files()
        out_files = []
        for f in self.files:
            if all([f[k]==v for k,v in kwa.iteritems()]):
                out_files.append(f)

        return out_files

    def init_gspread(self):
        self.gc = gspread.authorize(self.credentials)
