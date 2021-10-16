#usr/bin/env sh
python -c "import aiida; aiida.load_dbenv();from aiida.orm import WorkflowFactory;WorkChain = WorkflowFactory('elastic')"
