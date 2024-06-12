
How to add a new image
----------------------

Checkout ``pokeball/clean`` branch and navigate to ``azure-pipelines.yaml`` file. Add a new build step for your new image (see "AnalysisBuild" stage).

.. code-block:: yaml

   - stage: AnalysisBuild # please change this in your new stage
     displayName: 'Analysis Build' # also this
     variables: 
     - name: dockerfilePath
       value: '$(Build.SourcesDirectory)/Dockerfile_analysis_jupyter' # and this
     jobs:
     - template: 'ci-templates/imagebuild.yaml'
       parameters:
         jobName: 'BuildImage'
         dependsOn: ''
         condition: ''
         continueOnError: 'false'
         poolName: 'ADLS_OPENSTACK_POKEBALL'
         environment: 'test'
         folders: analysis # please change this too

To the deploy stage below add a set of new variables for instance "third_image", here is an example:

.. code-block:: yaml

    - name: third_image
       value: analysis:$(Build.SourceVersion)
     - name: third_imagename
       value: $(dockerRegistryServiceConnection).azurecr.io/jupyter/$(third_image)
     - name: third_imagetitle
       value: Pokeball-NewAmazingImage-latest

Please always count up, so if "third\ *..." is already taken, go ahead with "fourth*..." and so on.

After that, open ``ci-templates/deploy.yaml``\ and add your variables in the curl (currently line 21). Like this (see third_imagetitle, third_imagename):

.. code-block:: yaml

   ... "variables":{"first_imagetitle":{"value":"$(first_imagetitle)"},"first_imagename":{"value":"$(first_imagename)"},"second_imagetitle":{"value":"$(second_imagetitle)"},"third_imagetitle":{"value":"$(third_imagetitle)", "third_imagename":{"value":"$(third_imagename)"}}}'

This triggers the remote CaaS pipeline and your image will be then visable in the dropdown on cluster 6.

ATTENTION!! if you want to add a **fourth image or more** you need to currently still inform us (on Slack) so that new placeholder variables can be added on our side too.
Otherwise there would be a lot of ``placeholder`` in the dropdown. Working on a way to maybe make that still more user friendly.

----

Just for reference for us:

In the ``helm-charts/jupyterprofile/jupyterprofile.standalone.yml`` a new placeholder for each new image has to be added ("third_..." is already included).

These variables have to be added to the pipeline directly too:
Navigate to AOCC-Kubernetes-Remote, klick edit, then the button ``Variables`` at the top right and add the new variable pair there. "Let users override this value..." has to be selected and some placehoder value entered.

AOCC-Kubernetes-Remote works currently with the /templates/crd-update.yml file (TODO: I need to generalize it)
