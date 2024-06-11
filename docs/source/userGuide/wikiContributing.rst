
Wiki Contributing Guide
=======================

----

Thank you for investing your time in contributing to our `Kaiju project! <https://bitbucket.org/aplkaiju/kaiju>`_ :dragon_face:

Wiki Documentation
------------------

This wiki uses the `Markdown <http://daringfireball.net/projects/markdown/>`_ syntax. The `MarkDownDemo tutorial <https://bitbucket.org/tutorials/markdowndemo>`_ shows how various elements are rendered. The `Bitbucket documentation <https://confluence.atlassian.com/x/FA4zDQ>`_ has more information about using a wiki.

The wiki itself is actually a git repository, which means you can clone it, edit it locally/offline, add images or any other file type, and push it back to us. It will be live immediately.

Go ahead and try:

.. code-block::

   $ git clone https://bitbucket.org/aplkaiju/docs.git/wiki

Wiki pages are normal files, with the .md extension. You can edit them locally, as well as creating new ones.

**We have adopted the convention that all file and directory names will use `camelCase <https://en.wikipedia.org/wiki/Camel_case>`_.**

Syntax highlighting
-------------------

You can also highlight snippets of text (we use the excellent `Pygments <http://pygments.org/>`_ library).

Here's an example of some Python code:

.. code-block::

   #!python

   def wiki_rocks(text):
       formatter = lambda t: "funky"+t
       return formatter(text)

You can check out the source of this page to see how that's done, and make sure to bookmark `the vast library of Pygment lexers <http://pygments.org/docs/lexers/>`_\ , we accept the 'short name' or the 'mimetype' of anything in there.

Have fun!
