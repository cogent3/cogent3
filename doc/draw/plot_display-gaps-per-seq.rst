.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_draw_plot_display-gaps-per-seq.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_draw_plot_display-gaps-per-seq.py:


Counting gaps per sequence
==========================

We have several different ways of counting sequence gaps, and of visualising the results. By default, the `count_gaps_per_seq()` method returns a matrix of counts without the ability to visualise the results. When setting the argument `unique=True`, the counts are for gaps uniquely induced by each sequence. This can be a useful indicator of highly divergent sequences.


.. code-block:: default

    from cogent3 import load_aligned_seqs

    aln = load_aligned_seqs('../../tests/data/brca1.fasta', moltype='dna')

    counts = aln.count_gaps_per_seq(unique=True)
    counts






.. only:: builder_html

    .. raw:: html

        <table>
        <style>
        tr:last-child {border-bottom: 1px solid #000;} tr > th {text-align: center !important;} tr > td {text-align: left !important;}
        </style>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>FlyingFox</th>
        <th>DogFaced</th>
        <th>FreeTaile</th>
        <th>LittleBro</th>
        <th>TombBat</th>
        <th>RoundEare</th>
        <th>FalseVamp</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>LeafNose</th>
        <th>Horse</th>
        <th>Rhino</th>
        <th>Pangolin</th>
        <th>Cat</th>
        <th>Dog</th>
        <th>Llama</th>
        <th>Pig</th>
        <th>Cow</th>
        <th>Hippo</th>
        <th>SpermWhale</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">3</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>HumpbackW</th>
        <th>Mole</th>
        <th>Hedgehog</th>
        <th>TreeShrew</th>
        <th>FlyingLem</th>
        <th>Galago</th>
        <th>HowlerMon</th>
        <th>Rhesus</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">3</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">3</td>
        <td style="font-family: monospace !important;">21</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>Orangutan</th>
        <th>Gorilla</th>
        <th>Human</th>
        <th>Chimpanzee</th>
        <th>Jackrabbit</th>
        <th>FlyingSqu</th>
        <th>OldWorld</th>
        <th>Mouse</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">57</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>Rat</th>
        <th>NineBande</th>
        <th>HairyArma</th>
        <th>Anteater</th>
        <th>Sloth</th>
        <th>Dugong</th>
        <th>Manatee</th>
        <th>AfricanEl</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>AsianElep</th>
        <th>RockHyrax</th>
        <th>TreeHyrax</th>
        <th>Aardvark</th>
        <th>GoldenMol</th>
        <th>Madagascar</th>
        <th>Tenrec</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">6</td>
        </tr>
        </tbody>
        </table>
        <table>
        <thead style="background: rgba(161, 195, 209, 0.75); font-weight: bold; text-align: center;">
        <th>LesserEle</th>
        <th>GiantElep</th>
        <th>Caenolest</th>
        <th>Phascogale</th>
        <th>Wombat</th>
        <th>Bandicoot</th>
        </thead>
        <tbody>
        <tr>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">6</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        <td style="font-family: monospace !important;">0</td>
        </tr>
        </tbody>
        </table>

        <br />
        <br />

Plotting counts of unique gaps
##############################

There are three plot types supported. In all cases, placing the mouse pointer over a data point will show hover text with the sequence name.

Displaying unique gaps as a bar chart
*************************************


.. code-block:: default


    counts = aln.count_gaps_per_seq(unique=True, drawable='bar')
    counts.show(renderer="sphinx_gallery", width=500)



.. raw:: html
    :file: images/sphx_glr_plot_display-gaps-per-seq_001.html





Displaying unique gaps as a violin plot
***************************************


.. code-block:: default


    counts = aln.count_gaps_per_seq(unique=True, drawable='violin')
    counts.show(renderer="sphinx_gallery", width=300, height=500)



.. raw:: html
    :file: images/sphx_glr_plot_display-gaps-per-seq_002.html





Displaying unique gaps as a box plot
************************************


.. code-block:: default


    counts = aln.count_gaps_per_seq(unique=True, drawable='box')
    counts.show(renderer="sphinx_gallery", width=300, height=500)




.. raw:: html
    :file: images/sphx_glr_plot_display-gaps-per-seq_003.html






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  3.570 seconds)


.. _sphx_glr_download_draw_plot_display-gaps-per-seq.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: plot_display-gaps-per-seq.py <plot_display-gaps-per-seq.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: plot_display-gaps-per-seq.ipynb <plot_display-gaps-per-seq.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
