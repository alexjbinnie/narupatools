import testmodule

from apiparser.parse import parse_module

analyzed = parse_module(testmodule)

print(analyzed)