/* * GFLIP - Geometrical FLIRT Phrases for Large Scale Place Recognition
 * Copyright (C) 2012-2013 Gian Diego Tipaldi and Luciano Spinello and Wolfram
 * Burgard
 *
 * This file is part of GFLIP.
 *
 * GFLIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GFLIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GFLIP.  If not, see <http://www.gnu.org/licenses/>.
 */

/*! \mainpage 
 *
 * \section need Why do I need Geometrical FLIRT Phrases (GFLIP)?
 * Let's say you have a large 2D scan dataset and you need to localize, or to check for loop closing.
 * Usually, a brute forcing ICP scan matching one vs all solves the problem. 
 * In case you need to do this process online, or repetitively, this may be computationally expensive. 
 * 
 * GFLIP comes to help in this situations. GFLIP makes use of a FLIRT bag-of-word representation, to create "phrases": combinations of words that 
 * are robust to noise and are able to encode complicate sequential patterns. In contrast to other approaches, GFLIP exploits the fact that a scan
 *  is a 1D manifold to create a rotation invariant scan representation.
 * This allows for reliable and fast retreival of 2D laser scans. Way more reliable 
 * than standard bag-of-words approaches.
 * 
 * In practice, given a scan represented by FLIRT words eg:
 * 
 * \verbatim scan example: 50 12 56 87 94 52 68 \endverbatim
 * 
 * GFLIP matches it with a dataset containing 
 * \verbatim scan1: 13 12 56 8 94 52 687
scan2: 4 58 66 45 33 6 1
scan3: 1 223 3
scan4: 45 33 6 1 12 56
...  \endverbatim
 * 
 * \section intro Why are GFLIPs useful?
 *
 * Place recognition, i.e., the problem of recognizing
 * if the robot is navigating in an already visited place, is a
 * fundamental problem in mobile robot navigation. Efﬁcient
 * solutions to this problem are relevant for effectively localizing
 * robots and for creating maps in real time. Relatively few methods
 * have been proposed to efﬁciently solve this problem in very large
 * environments using 2D range data. <b>In this paper, we introduce
 * geometrical FLIRT phrases (GFLIPs) as a novel retrieval method
 * for very efﬁcient and precise place recognition. GFLIPs perform
 * approximate 2D range data matching, have low computational
 * cost, can handle complicated partial matching patterns and are
 * robust to noise</b>. Experiments carried out with publicly available
 * datasets demonstrate that GFLIPs largely outperform state-of-theart approaches in 2D range-based place recognition in terms of
 * efﬁciency and recall. We obtain retrieval performances with more
 * than 85% recall at 99% precision in less than a second, even on
 * data sets obtained from several kilometer long runs.
 * 
 * see the <a href="http://www.informatik.uni-freiburg.de/~spinello/tipaldiICRA13.pdf">paper</a> 
 * \section install_sec Installation
 * Download the library. The library relies on cmake to generate the Makefiles.
 * 
 * Go to the \c root directory of your project and run
 * \verbatim $ mkdir build
$ cd build  
$ cmake ../ 
$ make \endverbatim
 Binaries are generated in \c build/bin and \c build/lib
 * 
 * Library can be installed in your system by running
 * \verbatim $ make install \endverbatim
 * 
 * 
 * 
 * The software depends on the following external libraries
 * \li <em> Boost >= 1.4 (special_functions/binomial, serialization) </em>
 * \li <em> Eigen >= 3  </em>
 * \li <em> FLIRT (The cmake scripts take care of downloading it for you)  </em>
 *
 * \section ex_sec Binaries
 * Binaries can be found in \c bin/. For testing a FLIRT words dataset and the vocabulary used for the paper have been included in \c data/
 * The datasets can be found in hte data directory of the FLIRT library.
 * 
 * \verbatim gflip_cl \endverbatim reads a dataset where each scan is composed of FLIRT words and retreives best matches in the dataset for benchmarking. 
 * Many parameters can be selected from command line (kernel size, bag-of-words/distances, etc)
 *   GFPLoopClosingTest  learnVocabularyKMeans  
 *
 * \verbatim gflip_cl_onequery \endverbatim same as above but with one scan. This is an example for understanding how to use it.
 *
 * \verbatim learnVocabularyKMeans \endverbatim learns a vocabulary from a set of datasets.
 * 
 * \verbatim featureExtractor \endverbatim extracts the features from a dataset and writes them to a file for further use.
 * 
 * \verbatim generateBoW \endverbatim generates a BoW description of the whole dataset to be used for gflip_cl and gflip_cl_onequery.
 * 
 * \verbatim generateNN \endverbatim generates a nearest neighbor file to be used with nnLoopClosingTest.
 * 
 * \verbatim GFPLoopClosingTest \endverbatim performs the full matching from feature extraction to BoW generation and kernelized matching.
 * 
 * \section ref How to cite GFLIP
 *  "Geometrical FLIRT Phrases for Large Scale Place Recognition in 2D Range Data", G. D. Tipaldi, L. Spinello, W. Burgard -- Int. Conf. Robotics and Automation (ICRA) 2013
 *  
 * \section license License
 * GFLIP - Geometrical FLIRT Phrases for Large Scale Place Recognition
 * Copyright (C) 2012-2013 Gian Diego Tipaldi and Luciano Spinello and Wolfram Burgard
 *
 * This file is part of GFLIP.
 *
 * GFLIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GFLIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GFLIP.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */
 
