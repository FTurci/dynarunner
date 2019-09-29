import dynarunner


dynamo = dynarunner.DynamoMD("/Users/francesco/Downloads/dynamomd-1.7.128-OSX-amd64-clang/bin/")

dynamo.configure(0.5,20, polydispersity=0.08, N=6000)
dynarunner.conf_to_atom("output/config.start.xml",0)
dynamo.equilibrate(600*1000,conf=dynamo.start)
dynamo.compress(target_packing=0.46,conf=dynamo.equilibrated)
dynamo.equilibrate(600*100)
dynamo.run(dynamo.equilibrated, 100,10)
# dynarunner.conf_to_xyz("output/config.start.xml")
# dynarunner.conf_to_atom("output/config.start.xml",0)