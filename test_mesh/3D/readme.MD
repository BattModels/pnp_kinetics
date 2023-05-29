# Testing mesh for solution convergence

1. Set location of mesh (`.xml`) in `config.yml`.

2. Make bash files executable by running `chmod +x run.sh solve_vmg.sh process_vmg.sh`
3. Run `./run.sh Sol` : Runs PNP equation solve on a range of $`V_{app}`$ specified in `solve_vmg.sh`. Final currents and overpotentials are collected in `Sol/iv_data.txt`

2. `gmsh2xml.py` : Convert `.msh` output to dolfin-xml `.xml` readable by fenics 2019. 