import docker

client = docker.from_env()
image = 'drum:latest'
entrypoint = '/sw/drum/bin/trips-drum -nouser'
detach = True
internal_port = 6200
for port_inc in range(50):
    expose_port = 6201 + port_inc
    ports = {('%d/tcp' % internal_port): expose_port}
    cont = client.containers.run(image, entrypoint=entrypoint,
                                 detach=detach, ports=ports)
