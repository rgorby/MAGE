JupyterLab on NASA HEC Systems
==============================

What follows below is to compile and use Jupyterlab on Pleiades. Replace
rmalbarr with your specific Pleiades username. Note that you will have a
certain /nobackup directory number (for me it is /nobackupp12/rmalbarr).
Change this according to your number. Similarly, I have /home7/rmalbarr as my
home directory. Change this home directory number for you, accordingly.

Get set up on Pleiades
----------------------


Setup RSA tokens: copy from bitbucket to Pleiades Font End (pfe) and
`follow wiki <https://www.nas.nasa.gov/hecc/support/kb/enabling-your-rsa-securid-hard-token-(fob>`__59.html)

Setup of ssh pass through,
`follow wiki <https://www.nas.nasa.gov/hecc/support/kb/setting-up-ssh-passthrough_232.html>`_

This should enable you to login to pfe with where the passcode is given by
SecureID mobile app (dual-factor authentication):

.. code-block:: bash

   ssh pfe


Setup the sup client
`from wiki <https://www.nas.nasa.gov/hecc/support/kb/using-the-secure-unattended-proxy-(sup>`__145.html)

This will enable you to send large files between remote and local servers with
``shiftc``

From local machine, for example, run:

.. code-block:: bash

   sup shiftc <files> rmalbar@pfe:/nobackupp12/rmalbarr


Clone repo to your nobackup or home directory on pfe. Here, for example, check
out the hidra branch of kaiju. From pfe prompt, run:

.. code-block:: bash

   install home-brew
   git lsf install
   git branch
   git checkout hidra

From
`kaiju wiki <https://bitbucket.org/aplkaiju/kaiju/wiki/quickStart/prerequisites>`_,
run:

.. code-block:: bash

   git clone https://YOUR_BITBUCKET_USERNAME@bitbucket.org/aplkaiju/kaiju.git
   export KAIJUHOME=$HOME/kaiju

where $HOME, for me, is set to /nobackupp12/rmalbarr. Change this according to
your username and desired directory. Note: Use Atlassian App password for the
clone command above. This is found in your bitbucket profile.

Using Jupyter notebook on Pleiades
----------------------------------


Setup kaiju-3.8 virtual conda environment and Jupyter notebook kernel. From
`kaiju wiki <https://bitbucket.org/aplkaiju/kaiju/wiki/quickStart/install_python.md>`_
copy the .yml file called ``kaiju-3.8_pleiades_20220817.yml`` to your pfe home
directory. Edit the last line in the .yml file to your correct path:
``/swbuild/analytix/tools/miniconda3_220407/envs/kaiju-3.8``

From pfe home directory prompt, run:

.. code-block:: bash

   module use -a /swbuild/analytix/tools/modulefiles
   module load miniconda3/v4
   conda env create -f kaiju-3.8_pleiades_20220817.yml
   source activate kaiju-3.8
   pip install ipykernel
   python -m ipykernel install --user --name=kaiju-npl

If done correctly, it should output:

.. code-block:: bash

   Installed kernelspec kaiju-npl in /home7/rmalbarr/.local/share/jupyter/kernels/kaiju-npl

Note: You can check the available environments by running:

.. code-block:: bash

   conda info --envs


Reserve dedicated node for Jupyter analysis, where, below uses ivy system
(model=ivy) which is common for dedicated nodes. See
`wiki here for more information <https://www.nas.nasa.gov/hecc/support/kb/using-jupyter-notebook-for-machine-learning-development-on-nas-systems_576.html>`_
From pfe run an interactive mode and note the number of your pfe at command
prompt (which will change each pfe login), e.g. for me it is rmalbarr@pfe22.
From pfe run:

.. code-block:: bash

   qsub -I -lselect=1:ncpus=20:model=ivy,walltime=2:00:00

Once on a dedicated node, e.g., node number r437i4n1, navigate to your
directory that has the Jupyter notebook files, then run:

.. code-block:: bash

   ssh r437i4n1
   module use -a /swbuild/analytix/tools/modulefiles
   module load miniconda3/v4
   source activate jupyterlab
   jupyter lab --no-browser

From you local machine terminal, with the pfe number (e.g., here I have pfe22)
used to reserved the dedicated node and node number used to activate jupyter
lab (e.g., here I have r437i4n1), I run (for my username rmalbarr):

.. code-block:: bash

   ssh -l rmalbarr -o “StrictHostKeyChecking ask” -L 8080:localhost:8888 -o ProxyJump=sfe,pfe22 r437i4n1

Then, from your local machine browser, navigate to: https://localhost:8080/ or
https://127.0.0.1:8080 Here, you should be able to access the Jupyter notebook
files and the kaiju-3.8 kernel should be available and working.

Please reach out with any comments or questions to Robert Albarran, Ph.D. via
email to albarran1@atmos.ucla.edu
