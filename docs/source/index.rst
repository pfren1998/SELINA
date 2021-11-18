.. Selina documentation master file, created by
   sphinx-quickstart on Sat Nov  6 10:34:40 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Selina's documentation
=================================
Selina is a deep-learning based framework for single cell assignment with multiple references. The algorithm consists of three main steps: pre-training, fine-tuning and predicting. The reference datasets were first trained on MADA, a supervised deep learning framework, to obtain a pre-trained model. An autoencoder was then used to fine-tune the parameters of the pre-trained model. Finally, the labels from reference datasets were transferred to the query dataset based on the fully trained model. Along with the annotation algorithm, we also collected xx datasets which were uniformly processed and curated to provide users with comprehensive pre-trained models.

.. image:: _images/algorithm.png
   :width: 1600


Usage
=====


.. toctree::
   :maxdepth: 4

   prepare

.. toctree::
   :maxdepth: 4

   run
   



