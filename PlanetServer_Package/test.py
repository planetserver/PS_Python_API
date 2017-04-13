from PlanetServer import api

r = "R1330"
g = "D2300"
b = "HCPINDEX3"

# Note : Please take care that the coverage id has letters in lower case!
#coverage_id = "frt000064d9_07_if166l_trr3"
coverage_id = "frt000064d9_07_if166l_trr3"

print(api.coordinates_trans(50,10))
api.URL_creator(coverage_id, r, g, b)

