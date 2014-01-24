
! Focussing routine created by EXTRACTPKCN1.FOR
! Ian Harris 1994 (ibh11@cus.cam.ac.uk)
!
! GENIE routine to focus spectra from TEST Beam
!
! Unsmoothed background is not subtracted
! Final spectrum normalised by unsmoothed
!   Vanadium spectrum (#5138)
! Focussed to Detector 205
! Normalised to beam monitor count
!
! Template from KSK and JSW
!
$ nwork=p1
$ nrun=p2
$ nbin=p3
!
> ass 'nrun'
!
> set par w15 11.1116 0 180.0 0 0
> w15=s205
> u/d w15
> w'nwork'=w15
!
> set par w15 11.8733 0 180.0 0 0
> w15=s138
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8314 0 180.0 0 0
> w15=s139
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8015 0 180.0 0 0
> w15=s140
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7651 0 180.0 0 0
> w15=s141
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7356 0 180.0 0 0
> w15=s142
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6971 0 180.0 0 0
> w15=s143
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6659 0 180.0 0 0
> w15=s144
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6333 0 180.0 0 0
> w15=s145
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6029 0 180.0 0 0
> w15=s146
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5665 0 180.0 0 0
> w15=s147
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5314 0 180.0 0 0
> w15=s148
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4986 0 180.0 0 0
> w15=s149
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4656 0 180.0 0 0
> w15=s150
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4249 0 180.0 0 0
> w15=s151
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3938 0 180.0 0 0
> w15=s152
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3515 0 180.0 0 0
> w15=s153
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3217 0 180.0 0 0
> w15=s154
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2876 0 180.0 0 0
> w15=s155
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2547 0 180.0 0 0
> w15=s156
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2186 0 180.0 0 0
> w15=s157
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1837 0 180.0 0 0
> w15=s158
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1452 0 180.0 0 0
> w15=s159
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1133 0 180.0 0 0
> w15=s160
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0772 0 180.0 0 0
> w15=s161
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0410 0 180.0 0 0
> w15=s162
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0000 0 180.0 0 0
> w15=s163
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9710 0 180.0 0 0
> w15=s164
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9253 0 180.0 0 0
> w15=s165
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8995 0 180.0 0 0
> w15=s166
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8627 0 180.0 0 0
> w15=s167
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8253 0 180.0 0 0
> w15=s168
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7814 0 180.0 0 0
> w15=s169
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7501 0 180.0 0 0
> w15=s170
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7099 0 180.0 0 0
> w15=s171
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6779 0 180.0 0 0
> w15=s172
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6388 0 180.0 0 0
> w15=s173
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6007 0 180.0 0 0
> w15=s174
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5597 0 180.0 0 0
> w15=s175
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5282 0 180.0 0 0
> w15=s176
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4859 0 180.0 0 0
> w15=s177
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4532 0 180.0 0 0
> w15=s178
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4133 0 180.0 0 0
> w15=s179
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3803 0 180.0 0 0
> w15=s180
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3346 0 180.0 0 0
> w15=s181
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8800 0 180.0 0 0
> w15=s183
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8404 0 180.0 0 0
> w15=s184
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8069 0 180.0 0 0
> w15=s185
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7720 0 180.0 0 0
> w15=s186
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7415 0 180.0 0 0
> w15=s187
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7018 0 180.0 0 0
> w15=s188
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6731 0 180.0 0 0
> w15=s189
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6349 0 180.0 0 0
> w15=s190
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6065 0 180.0 0 0
> w15=s191
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5690 0 180.0 0 0
> w15=s192
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5374 0 180.0 0 0
> w15=s193
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4971 0 180.0 0 0
> w15=s194
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4664 0 180.0 0 0
> w15=s195
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4247 0 180.0 0 0
> w15=s196
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3944 0 180.0 0 0
> w15=s197
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3568 0 180.0 0 0
> w15=s198
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3233 0 180.0 0 0
> w15=s199
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2868 0 180.0 0 0
> w15=s200
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2547 0 180.0 0 0
> w15=s201
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2134 0 180.0 0 0
> w15=s202
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1852 0 180.0 0 0
> w15=s203
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1460 0 180.0 0 0
> w15=s204
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0746 0 180.0 0 0
> w15=s206
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0404 0 180.0 0 0
> w15=s207
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9997 0 180.0 0 0
> w15=s208
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9688 0 180.0 0 0
> w15=s209
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9269 0 180.0 0 0
> w15=s210
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8960 0 180.0 0 0
> w15=s211
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8573 0 180.0 0 0
> w15=s212
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8168 0 180.0 0 0
> w15=s213
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7788 0 180.0 0 0
> w15=s214
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7452 0 180.0 0 0
> w15=s215
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7031 0 180.0 0 0
> w15=s216
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6743 0 180.0 0 0
> w15=s217
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6391 0 180.0 0 0
> w15=s218
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6021 0 180.0 0 0
> w15=s219
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5640 0 180.0 0 0
> w15=s220
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5232 0 180.0 0 0
> w15=s221
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4777 0 180.0 0 0
> w15=s222
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4506 0 180.0 0 0
> w15=s223
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4137 0 180.0 0 0
> w15=s224
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3725 0 180.0 0 0
> w15=s225
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3258 0 180.0 0 0
> w15=s226
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.2832 0 180.0 0 0
> w15=s227
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.8471 0 180.0 0 0
> w15=s229
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 85.9232 0 180.0 0 0
> w15=s230
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.7774 0 180.0 0 0
> w15=s231
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 ******* 0 180.0 0 0
> w15=s232
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6993 0 180.0 0 0
> w15=s233
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6660 0 180.0 0 0
> w15=s234
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6343 0 180.0 0 0
> w15=s235
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.6060 0 180.0 0 0
> w15=s236
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5732 0 180.0 0 0
> w15=s237
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.5391 0 180.0 0 0
> w15=s238
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4951 0 180.0 0 0
> w15=s239
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4700 0 180.0 0 0
> w15=s240
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.4271 0 180.0 0 0
> w15=s241
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3999 0 180.0 0 0
> w15=s242
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3569 0 180.0 0 0
> w15=s243
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.3320 0 180.0 0 0
> w15=s244
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2966 0 180.0 0 0
> w15=s245
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2602 0 180.0 0 0
> w15=s246
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.2233 0 180.0 0 0
> w15=s247
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1899 0 180.0 0 0
> w15=s248
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1398 0 180.0 0 0
> w15=s249
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.1124 0 180.0 0 0
> w15=s250
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0840 0 180.0 0 0
> w15=s251
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0509 0 180.0 0 0
> w15=s252
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 11.0090 0 180.0 0 0
> w15=s253
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9762 0 180.0 0 0
> w15=s254
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9346 0 180.0 0 0
> w15=s255
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.9025 0 180.0 0 0
> w15=s256
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8659 0 180.0 0 0
> w15=s257
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.8327 0 180.0 0 0
> w15=s258
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7924 0 180.0 0 0
> w15=s259
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.7221 0 180.0 0 0
> w15=s261
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6852 0 180.0 0 0
> w15=s262
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6484 0 180.0 0 0
> w15=s263
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.6115 0 180.0 0 0
> w15=s264
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5700 0 180.0 0 0
> w15=s265
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.5349 0 180.0 0 0
> w15=s266
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4949 0 180.0 0 0
> w15=s267
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4667 0 180.0 0 0
> w15=s268
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.4248 0 180.0 0 0
> w15=s269
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3858 0 180.0 0 0
> w15=s270
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.3409 0 180.0 0 0
> w15=s271
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
> set par w15 10.2980 0 180.0 0 0
> w15=s272
> u/d w15
> rebin w'nwork' w15
> w'nwork'=w'nwork'+w15
!
! ***** NORMALIZE TO INCIDENT BEAM MONITOR *****
!
> w15=s1
> integrate w15 5000 5400
> W'NWORK'=W'NWORK'/V1
!
> SET PAR W'NWORK' 15.0 0.7142 90 0 0
> SET YCODE W'NWORK' -1
> U/T W'NWORK'
!
! ***** SUBTRACT INSTRUMENTAL BACKGROUND *****
!
!> read/open w15 teb_calib_empty:ef5137.dat
!> read/close
!> rebin w'nwork' w15
!> w'nwork'=w'nwork'-w15
!
! ***** NORMALIZE TO VANADIUM SPECTRUM *****
!
> read/open w15 teb_calib_van:ef223054.his 2
> read/close
> rebin w'nwork' W15
> w'nwork'=w'nwork'/w15
!
> rebin w'nwork' 1000:19000
!
> a b 'nbin'
> u/t w'nwork'
!