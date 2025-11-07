import jawm

featurecounts_submodule_p1=jawm.Process( 
    name="featurecounts_submodule_p1",
    script="""#!/bin/bash
echo "Demo module echo process"
"""  
)