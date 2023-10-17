# PreComp

Classical laminate theory for composite wings/blades based on: https://nwtc.nrel.gov/PreComp

Please make all feature changes and bug fixes as branches and then create pull requests against the dev branch.  The dev branch will be periodically pulled into master for version changes.

This code has been translated to the Julia programming language and modified to handle automatic differentiation by Taylor McDonnell at Brigham Young University during his PhD work.

Even though this code has been publicly avalaible on the NREL website in the past, is extensively used in the wind energy sector, and is currently available at https://wind.nrel.gov/forum/wind/viewtopic.php?t=2575, from email communication with Jason Jonkman at NREL, apparently it never went through a formal release process at NREL but is covered by their disclaimer https://www.nrel.gov/wind/nwtc/disclaimer.html.  The original code reserves all rights against modification or use of the code in its header - "This code is the property of NREL.  The code, or any part of it, may not be transferred, altered or copied, directly or indirectly, without approval from NWTC/NREL."

Since the code has been distributed in the public domain for quite some time (years), it appears reasonable to assume that NREL has granted the approval to transfer, alter, copy, directly and indirectly, the code, especially considering the support tickets linked above, and the code's extensive use in the public domain.  However, all credit for the translation and improvements of this version should be given to Taylor McDonnell.  It should also be noted that the original disclaimer still applies.
