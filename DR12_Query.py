"""
from astroquery.SciServer import SciServer
from SciServer import CasJobs
from SciServer import LoginPortal
# Read the current access token from the file

my_access_token = LoginPortal.getToken()
query = "SELECT TOP 10 p.objid, p.ra, p.dec, s.z, s.zerr, pz.z, pz.zerr FROM PhotoObj p, SpecObjAll s, Photoz pz WHERE p.dec >=-53 AND p.dec<= 90 "
responseStream = CasJobs.executeQuery(query, "DR12", token =
        my_access_token)
result = responseStream.read()

print("\n---Query---\n{}\n---Result---\n{}".format(query, result))
