#define I_MOY 0
#define UI_MOY 1
#define VI_MOY 2
#define WI_MOY 3
#define PI_MOY 4
#define UIv_MOY 5
#define VIv_MOY 6
#define WIv_MOY 7
#define PIv_MOY 8
// Correlations vitesse :
// 2eme ordre :
#define UUI_MOY 9
#define VVI_MOY 10
#define WWI_MOY 11
#define UVI_MOY 12
#define VWI_MOY 13
#define UWI_MOY 14
// 3eme ordre :
#define UUUI_MOY 15
#define UUVI_MOY 16
#define UUWI_MOY 17
#define UVVI_MOY 18
#define UVWI_MOY 19
#define UWWI_MOY 20
#define VVVI_MOY 21
#define VVWI_MOY 22
#define VWWI_MOY 23
#define WWWI_MOY 24
//
// 2eme ordre vit/pression :
#define UPI_MOY 25
#define VPI_MOY 26
#define WPI_MOY 27
//
//
// Maintenant, ce qui fait intervenir les derivees :
// 1. de l'indicatrice  :
//^^^^^^^^^^^^^^^^^^^^^^^
// ui*dI/dxb : (9)
#define UdIdx_MOY 28
#define VdIdx_MOY 29
#define WdIdx_MOY 30
#define UdIdy_MOY 31
#define VdIdy_MOY 32
#define WdIdy_MOY 33
#define UdIdz_MOY 34
#define VdIdz_MOY 35
#define WdIdz_MOY 36
//^^^^^^^^^^^^^^^^^^^^^^^
// ui*p*dI/dxb : (9)
#define UPdIdx_MOY 37
#define VPdIdx_MOY 38
#define WPdIdx_MOY 39
#define UPdIdy_MOY 40
#define VPdIdy_MOY 41
#define WPdIdy_MOY 42
#define UPdIdz_MOY 43
#define VPdIdz_MOY 44
#define WPdIdz_MOY 45
//^^^^^^^^^^^^^^^^^^^^^^^
// ui*uj*dI/dxb : (18)
#define UUdIdx_MOY 46
#define VVdIdx_MOY 47
#define WWdIdx_MOY 48
#define UUdIdy_MOY 49
#define VVdIdy_MOY 50
#define WWdIdy_MOY 51
#define UUdIdz_MOY 52
#define VVdIdz_MOY 53
#define WWdIdz_MOY 54
#define UVdIdx_MOY 55
#define UWdIdx_MOY 56
#define VWdIdx_MOY 57
#define UVdIdy_MOY 58
#define UWdIdy_MOY 59
#define VWdIdy_MOY 60
#define UVdIdz_MOY 61
#define UWdIdz_MOY 62
#define VWdIdz_MOY 63
//^^^^^^^^^^^^^^^^^^^^^^
// p*dI/dxb : (3)
#define PdIdx_MOY 64
#define PdIdy_MOY 65
#define PdIdz_MOY 66
//^^^^^^^^^^^^^^^^^^^^^^
// I*dI/dxb : (3)
#define IdIdx_MOY 67
#define IdIdy_MOY 68
#define IdIdz_MOY 69
// 2. de la vitesse  :
//^^^^^^^^^^^^^^^^^^^^^^^
// I*dUidxb : (9)
#define IdUdx_MOY 70
#define IdVdx_MOY 71
#define IdWdx_MOY 72
#define IdUdy_MOY 73
#define IdVdy_MOY 74
#define IdWdy_MOY 75
#define IdUdz_MOY 76
#define IdVdz_MOY 77
#define IdWdz_MOY 78
// I*P*dUidxb : (9)
#define IPdUdx_MOY 79
#define IPdVdx_MOY 80
#define IPdWdx_MOY 81
#define IPdUdy_MOY 82
#define IPdVdy_MOY 83
#define IPdWdy_MOY 84
#define IPdUdz_MOY 85
#define IPdVdz_MOY 86
#define IPdWdz_MOY 87
// 3. de l'indicatrice et de la vitesse  :
//^^^^^^^^^^^^^^^^^^^^^^^
// dUidxb*dIdxb : (9)
#define dUdxdIdx_MOY 88
#define dVdxdIdx_MOY 89
#define dWdxdIdx_MOY 90
#define dUdydIdy_MOY 91
#define dVdydIdy_MOY 92
#define dWdydIdy_MOY 93
#define dUdzdIdz_MOY 94
#define dVdzdIdz_MOY 95
#define dWdzdIdz_MOY 96
// Ui*dUjdxb*dIdxb : (27)
#define UdUdxdIdx_MOY 97
#define UdVdxdIdx_MOY 98
#define UdWdxdIdx_MOY 99
#define UdUdydIdy_MOY 100
#define UdVdydIdy_MOY 101
#define UdWdydIdy_MOY 102
#define UdUdzdIdz_MOY 103
#define UdVdzdIdz_MOY 104
#define UdWdzdIdz_MOY 105
//
#define VdUdxdIdx_MOY 106
#define VdVdxdIdx_MOY 107
#define VdWdxdIdx_MOY 108
#define VdUdydIdy_MOY 109
#define VdVdydIdy_MOY 110
#define VdWdydIdy_MOY 111
#define VdUdzdIdz_MOY 112
#define VdVdzdIdz_MOY 113
#define VdWdzdIdz_MOY 114
//
#define WdUdxdIdx_MOY 115
#define WdVdxdIdx_MOY 116
#define WdWdxdIdx_MOY 117
#define WdUdydIdy_MOY 118
#define WdVdydIdy_MOY 119
#define WdWdydIdy_MOY 120
#define WdUdzdIdz_MOY 121
#define WdVdzdIdz_MOY 122
#define WdWdzdIdz_MOY 123
//
// 4. de la vitesse et de la vitesse (fois l'indic) :
//^^^^^^^^^^^^^^^^^^^^^^^
// dUidxj*dUmdxb : (36+9 = 45)
#define IdUdxdUdx_MOY 124
#define IdUdxdUdy_MOY 125
#define IdUdxdUdz_MOY 126
#define IdUdxdVdx_MOY 127
#define IdUdxdVdy_MOY 128
#define IdUdxdVdz_MOY 129
#define IdUdxdWdx_MOY 130
#define IdUdxdWdy_MOY 131
#define IdUdxdWdz_MOY 132
#define IdUdydUdy_MOY 133
#define IdUdydUdz_MOY 134
#define IdUdydVdx_MOY 135
#define IdUdydVdy_MOY 136
#define IdUdydVdz_MOY 137
#define IdUdydWdx_MOY 138
#define IdUdydWdy_MOY 139
#define IdUdydWdz_MOY 140
#define IdUdzdUdz_MOY 141
#define IdUdzdVdx_MOY 142
#define IdUdzdVdy_MOY 143
#define IdUdzdVdz_MOY 144
#define IdUdzdWdx_MOY 145
#define IdUdzdWdy_MOY 146
#define IdUdzdWdz_MOY 147
#define IdVdxdVdx_MOY 148
#define IdVdxdVdy_MOY 149
#define IdVdxdVdz_MOY 150
#define IdVdxdWdx_MOY 151
#define IdVdxdWdy_MOY 152
#define IdVdxdWdz_MOY 153
#define IdVdydVdy_MOY 154
#define IdVdydVdz_MOY 155
#define IdVdydWdx_MOY 156
#define IdVdydWdy_MOY 157
#define IdVdydWdz_MOY 158
#define IdVdzdVdz_MOY 159
#define IdVdzdWdx_MOY 160
#define IdVdzdWdy_MOY 161
#define IdVdzdWdz_MOY 162
#define IdWdxdWdx_MOY 163
#define IdWdxdWdy_MOY 164
#define IdWdxdWdz_MOY 165
#define IdWdydWdy_MOY 166
#define IdWdydWdz_MOY 167
#define IdWdzdWdz_MOY 168
//
// 5. D'autre gradients plus regulier pour eviter d'utiliser le gradI :
//^^^^^^^^^^^^^^^^^^^^^^^
// I*dPdxb : (3)
#define IdPdx_MOY 169
#define IdPdy_MOY 170
#define IdPdz_MOY 171
//  I*Ui*dPdxb : (9)
#define IUdPdx_MOY 172
#define IUdPdy_MOY 173
#define IUdPdz_MOY 174
#define IVdPdx_MOY 175
#define IVdPdy_MOY 176
#define IVdPdz_MOY 177
#define IWdPdx_MOY 178
#define IWdPdy_MOY 179
#define IWdPdz_MOY 180
//  I*Ui*dUjdxb : (27)
#define IUdUdx_MOY 181
#define IUdVdx_MOY 182
#define IUdWdx_MOY 183
#define IUdUdy_MOY 184
#define IUdVdy_MOY 185
#define IUdWdy_MOY 186
#define IUdUdz_MOY 187
#define IUdVdz_MOY 188
#define IUdWdz_MOY 189
//
#define IVdUdx_MOY 190
#define IVdVdx_MOY 191
#define IVdWdx_MOY 192
#define IVdUdy_MOY 193
#define IVdVdy_MOY 194
#define IVdWdy_MOY 195
#define IVdUdz_MOY 196
#define IVdVdz_MOY 197
#define IVdWdz_MOY 198
//
#define IWdUdx_MOY 199
#define IWdVdx_MOY 200
#define IWdWdx_MOY 201
#define IWdUdy_MOY 202
#define IWdVdy_MOY 203
#define IWdWdy_MOY 204
#define IWdUdz_MOY 205
#define IWdVdz_MOY 206
#define IWdWdz_MOY 207
// I*ddUidxbxb : (9)
#define IddUdxx_MOY 208
#define IddVdxx_MOY 209
#define IddWdxx_MOY 210
#define IddUdyy_MOY 211
#define IddVdyy_MOY 212
#define IddWdyy_MOY 213
#define IddUdzz_MOY 214
#define IddVdzz_MOY 215
#define IddWdzz_MOY 216
//
//
// Other quantities :
//
// Force : pot * dIdxb
#define Fx_MOY 217
#define Fy_MOY 218
#define Fz_MOY 219
// Repulsion : pot_rep * dIdxb
#define Frx_MOY 220
#define Fry_MOY 221
#define Frz_MOY 222
//
// Repulsion : sum std::fabs(pot_rep * dIdxb)
#define Frax_MOY 223
#define Fray_MOY 224
#define Fraz_MOY 225
//
// frep frep_abs et f_total
//
// La dissipation :
#define DISSIP_MOY 226
//
// Second ordre vitesse (Cote vapeur)
// 2eme ordre :
#define UUIv_MOY 227
#define VVIv_MOY 228
#define WWIv_MOY 229
#define UVIv_MOY 230
#define VWIv_MOY 231
#define UWIv_MOY 232
// Autre interpolation :
#define UUIbv_MOY 233
#define VVIbv_MOY 234
#define WWIbv_MOY 235
#define UUIb_MOY 236
#define VVIb_MOY 237
#define WWIb_MOY 238
// 3eme ordre :
#define UUUIb_MOY 239
#define UUVIb_MOY 240
#define UUWIb_MOY 241
#define UVVIb_MOY 242
//  UVWIb_MOY n'existe pas...
#define UWWIb_MOY 243
#define VVVIb_MOY 244
#define VVWIb_MOY 245
#define VWWIb_MOY 246
#define WWWIb_MOY 247
//
// ai :
#define AI_MOY 248
#define NK_MOY 249
//
// Des correl supplementaires pour l'equation d'epsilon : (Are these the terms for the transport equation for the dissipation epsilon??
// ATTENTION! Leurs calculs ne sont pas tous codes!!
#define IdUdxdUdxdUdx_MOY 250
#define IdUdydUdydUdx_MOY 251
#define IdUdzdUdzdUdx_MOY 252
#define IdUdxdVdxdUdy_MOY 253
#define IdUdydVdydUdy_MOY 254
#define IdUdzdVdzdUdy_MOY 255
#define IdUdxdWdxdUdz_MOY 256
#define IdUdydWdydUdz_MOY 257
#define IdUdzdWdzdUdz_MOY 258
#define IdVdxdUdxdVdx_MOY 259
#define IdVdydUdydVdx_MOY 260
#define IdVdzdUdzdVdx_MOY 261
#define IdVdxdVdxdVdy_MOY 262
#define IdVdydVdydVdy_MOY 263
#define IdVdzdVdzdVdy_MOY 264
#define IdVdxdWdxdVdz_MOY 265
#define IdVdydWdydVdz_MOY 266
#define IdVdzdWdzdVdz_MOY 267
#define IdWdxdUdxdWdx_MOY 268
#define IdWdydUdydWdx_MOY 269
#define IdWdzdUdzdWdx_MOY 270
#define IdWdxdVdxdWdy_MOY 271
#define IdWdydVdydWdy_MOY 272
#define IdWdzdVdzdWdy_MOY 273
#define IdWdxdWdxdWdz_MOY 274
#define IdWdydWdydWdz_MOY 275
#define IdWdzdWdzdWdz_MOY 276
#define IdUdxdUdxW_MOY 277
#define IdUdydUdyW_MOY 278
#define IdUdzdUdzW_MOY 279
#define IdVdxdVdxW_MOY 280
#define IdVdydVdyW_MOY 281
#define IdVdzdVdzW_MOY 282
#define IdWdxdWdxW_MOY 283
#define IdWdydWdyW_MOY 284
#define IdWdzdWdzW_MOY 285
#define IdUdxddPdxdx_MOY 286
#define IdUdyddPdxdy_MOY 287
#define IdUdzddPdxdz_MOY 288
#define IdVdxddPdydx_MOY 289
#define IdVdyddPdydy_MOY 290
#define IdVdzddPdydz_MOY 291
#define IdWdxddPdzdx_MOY 292
#define IdWdyddPdzdy_MOY 293
#define IdWdzddPdzdz_MOY 294
#define IddUdxdxddUdxdx_MOY 295
#define IddUdxdyddUdxdy_MOY 296
#define IddUdxdzddUdxdz_MOY 297
#define IddUdydxddUdydx_MOY 298
#define IddUdydyddUdydy_MOY 299
#define IddUdydzddUdydz_MOY 300
#define IddUdzdxddUdzdx_MOY 301
#define IddUdzdyddUdzdy_MOY 302
#define IddUdzdzddUdzdz_MOY 303
#define IddVdxdxddVdxdx_MOY 304
#define IddVdxdyddVdxdy_MOY 305
#define IddVdxdzddVdxdz_MOY 306
#define IddVdydxddVdydx_MOY 307
#define IddVdydyddVdydy_MOY 308
#define IddVdydzddVdydz_MOY 309
#define IddVdzdxddVdzdx_MOY 310
#define IddVdzdyddVdzdy_MOY 311
#define IddVdzdzddVdzdz_MOY 312
#define IddWdxdxddWdxdx_MOY 313
#define IddWdxdyddWdxdy_MOY 314
#define IddWdxdzddWdxdz_MOY 315
#define IddWdydxddWdydx_MOY 316
#define IddWdydyddWdydy_MOY 317
#define IddWdydzddWdydz_MOY 318
#define IddWdzdxddWdzdx_MOY 319
#define IddWdzdyddWdzdy_MOY 320
#define IddWdzdzddWdzdz_MOY 321
#define dIdxddUdxdxdUdx_MOY 322
#define dIdxddUdxdydUdy_MOY 323
#define dIdxddUdxdzdUdz_MOY 324
#define dIdyddUdydxdUdx_MOY 325
#define dIdyddUdydydUdy_MOY 326
#define dIdyddUdydzdUdz_MOY 327
#define dIdzddUdzdxdUdx_MOY 328
#define dIdzddUdzdydUdy_MOY 329
#define dIdzddUdzdzdUdz_MOY 330
#define dIdxddVdxdxdVdx_MOY 331
#define dIdxddVdxdydVdy_MOY 332
#define dIdxddVdxdzdVdz_MOY 333
#define dIdyddVdydxdVdx_MOY 334
#define dIdyddVdydydVdy_MOY 335
#define dIdyddVdydzdVdz_MOY 336
#define dIdzddVdzdxdVdx_MOY 337
#define dIdzddVdzdydVdy_MOY 338
#define dIdzddVdzdzdVdz_MOY 339
#define dIdxddWdxdxdWdx_MOY 340
#define dIdxddWdxdydWdy_MOY 341
#define dIdxddWdxdzdWdz_MOY 342
#define dIdyddWdydxdWdx_MOY 343
#define dIdyddWdydydWdy_MOY 344
#define dIdyddWdydzdWdz_MOY 345
#define dIdzddWdzdxdWdx_MOY 346
#define dIdzddWdzdydWdy_MOY 347
#define dIdzddWdzdzdWdz_MOY 348
#define dIdxddUdxdz_MOY 349
#define dIdyddUdydz_MOY 350
#define dIdzddUdzdz_MOY 351
#define dIdzdUdxdUdx_MOY 352
#define dIdzdUdydUdy_MOY 353
#define dIdzdUdzdUdz_MOY 354
#define dIdzdVdxdVdx_MOY 355
#define dIdzdVdydVdy_MOY 356
#define dIdzdVdzdVdz_MOY 357
#define dIdzdWdxdWdx_MOY 358
#define dIdzdWdydWdy_MOY 359
#define dIdzdWdzdWdz_MOY 360
//
#define Nx_MOY 361
#define Ny_MOY 362
#define Nz_MOY 363
//
// Ajout de tout ce qui avait du dIdx...
// Mais retrait des doublons dans les derivees seconde :
// car ddUdydx = ddUdxdy (9 termes)
#define UaiNx_MOY 364
#define VaiNx_MOY 365
#define WaiNx_MOY 366
#define UaiNy_MOY 367
#define VaiNy_MOY 368
#define WaiNy_MOY 369
#define UaiNz_MOY 370
#define VaiNz_MOY 371
#define WaiNz_MOY 372
#define UPaiNx_MOY 373
#define VPaiNx_MOY 374
#define WPaiNx_MOY 375
#define UPaiNy_MOY 376
#define VPaiNy_MOY 377
#define WPaiNy_MOY 378
#define UPaiNz_MOY 379
#define VPaiNz_MOY 380
#define WPaiNz_MOY 381
#define UUaiNx_MOY 382
#define VVaiNx_MOY 383
#define WWaiNx_MOY 384
#define UUaiNy_MOY 385
#define VVaiNy_MOY 386
#define WWaiNy_MOY 387
#define UUaiNz_MOY 388
#define VVaiNz_MOY 389
#define WWaiNz_MOY 390
#define UVaiNx_MOY 391
#define UWaiNx_MOY 392
#define VWaiNx_MOY 393
#define UVaiNy_MOY 394
#define UWaiNy_MOY 395
#define VWaiNy_MOY 396
#define UVaiNz_MOY 397
#define UWaiNz_MOY 398
#define VWaiNz_MOY 399
#define PaiNx_MOY 400
#define PaiNy_MOY 401
#define PaiNz_MOY 402
#define IaiNx_MOY 403
#define IaiNy_MOY 404
#define IaiNz_MOY 405
#define dUdxaiNx_MOY 406
#define dVdxaiNx_MOY 407
#define dWdxaiNx_MOY 408
#define dUdyaiNy_MOY 409
#define dVdyaiNy_MOY 410
#define dWdyaiNy_MOY 411
#define dUdzaiNz_MOY 412
#define dVdzaiNz_MOY 413
#define dWdzaiNz_MOY 414
#define UdUdxaiNx_MOY 415
#define UdVdxaiNx_MOY 416
#define UdWdxaiNx_MOY 417
#define UdUdyaiNy_MOY 418
#define UdVdyaiNy_MOY 419
#define UdWdyaiNy_MOY 420
#define UdUdzaiNz_MOY 421
#define UdVdzaiNz_MOY 422
#define UdWdzaiNz_MOY 423
#define VdUdxaiNx_MOY 424
#define VdVdxaiNx_MOY 425
#define VdWdxaiNx_MOY 426
#define VdUdyaiNy_MOY 427
#define VdVdyaiNy_MOY 428
#define VdWdyaiNy_MOY 429
#define VdUdzaiNz_MOY 430
#define VdVdzaiNz_MOY 431
#define VdWdzaiNz_MOY 432
#define WdUdxaiNx_MOY 433
#define WdVdxaiNx_MOY 434
#define WdWdxaiNx_MOY 435
#define WdUdyaiNy_MOY 436
#define WdVdyaiNy_MOY 437
#define WdWdyaiNy_MOY 438
#define WdUdzaiNz_MOY 439
#define WdVdzaiNz_MOY 440
#define WdWdzaiNz_MOY 441
#define aiNxddUdxdxdUdx_MOY 442
#define aiNxddUdxdydUdy_MOY 443
#define aiNxddUdxdzdUdz_MOY 444
#define aiNyddUdydydUdy_MOY 445
#define aiNyddUdydzdUdz_MOY 446
#define aiNzddUdzdzdUdz_MOY 447
#define aiNxddVdxdxdVdx_MOY 448
#define aiNxddVdxdydVdy_MOY 449
#define aiNxddVdxdzdVdz_MOY 450
#define aiNyddVdydydVdy_MOY 451
#define aiNyddVdydzdVdz_MOY 452
#define aiNzddVdzdzdVdz_MOY 453
#define aiNxddWdxdxdWdx_MOY 454
#define aiNxddWdxdydWdy_MOY 455
#define aiNxddWdxdzdWdz_MOY 456
#define aiNyddWdydydWdy_MOY 457
#define aiNyddWdydzdWdz_MOY 458
#define aiNzddWdzdzdWdz_MOY 459
#define aiNxddUdxdz_MOY 460
#define aiNyddUdydz_MOY 461
#define aiNzddUdzdz_MOY 462
#define aiNzdUdxdUdx_MOY 463
#define aiNzdUdydUdy_MOY 464
#define aiNzdUdzdUdz_MOY 465
#define aiNzdVdxdVdx_MOY 466
#define aiNzdVdydVdy_MOY 467
#define aiNzdVdzdVdz_MOY 468
#define aiNzdWdxdWdx_MOY 469
#define aiNzdWdydWdy_MOY 470
#define aiNzdWdzdWdz_MOY 471
//
#define aiNx_MOY 472
#define aiNy_MOY 473
#define aiNz_MOY 474
#define kaiNx_MOY 475
#define kaiNy_MOY 476
#define kaiNz_MOY 477
//
#define aiNxx_MOY 478
#define aiNxy_MOY 479
#define aiNxz_MOY 480
#define aiNyy_MOY 481
#define aiNyz_MOY 482
#define aiNzz_MOY 483
#define aiNNxxxx_MOY 484
#define aiNNxxxy_MOY 485
#define aiNNxxxz_MOY 486
#define aiNNxxyy_MOY 487
#define aiNNxxyz_MOY 488
#define aiNNxxzz_MOY 489
#define aiNNxyyy_MOY 490
#define aiNNxyyz_MOY 491
#define aiNNxyzz_MOY 492
#define aiNNxzyy_MOY 493
#define aiNNxzyz_MOY 494
#define aiNNxzzz_MOY 495
#define aiNNyyyy_MOY 496
#define aiNNyyyz_MOY 497
#define aiNNyyzz_MOY 498
#define aiNNyzzz_MOY 499
#define aiNNzzzz_MOY 500
#define kai_MOY 501
#define UI_INT 502
#define VI_INT 503
#define WI_INT 504
#define UI_INTUI_INT 505
#define UI_INTVI_INT 506
#define UI_INTWI_INT 507
#define VI_INTVI_INT 508
#define VI_INTWI_INT 509
#define WI_INTWI_INT 510
#define UI_INTUI_INTUI_INT 511
#define UI_INTUI_INTVI_INT 512
#define UI_INTUI_INTWI_INT 513
#define UI_INTVI_INTUI_INT 514
#define UI_INTVI_INTVI_INT 515
#define UI_INTVI_INTWI_INT 516
#define UI_INTWI_INTUI_INT 517
#define UI_INTWI_INTVI_INT 518
#define UI_INTWI_INTWI_INT 519
#define VI_INTUI_INTUI_INT 520
#define VI_INTUI_INTVI_INT 521
#define VI_INTUI_INTWI_INT 522
#define VI_INTVI_INTUI_INT 523
#define VI_INTVI_INTVI_INT 524
#define VI_INTVI_INTWI_INT 525
#define VI_INTWI_INTUI_INT 526
#define VI_INTWI_INTVI_INT 527
#define VI_INTWI_INTWI_INT 528
#define WI_INTUI_INTUI_INT 529
#define WI_INTUI_INTVI_INT 530
#define WI_INTUI_INTWI_INT 531
#define WI_INTVI_INTUI_INT 532
#define WI_INTVI_INTVI_INT 533
#define WI_INTVI_INTWI_INT 534
#define WI_INTWI_INTUI_INT 535
#define WI_INTWI_INTVI_INT 536
#define WI_INTWI_INTWI_INT 537
#define UI_INTUIUI 538
#define UI_INTUIVI 539
#define UI_INTUIWI 540
#define UI_INTVIUI 541
#define UI_INTVIVI 542
#define UI_INTVIWI 543
#define UI_INTWIUI 544
#define UI_INTWIVI 545
#define UI_INTWIWI 546
#define VI_INTUIUI 547
#define VI_INTUIVI 548
#define VI_INTUIWI 549
#define VI_INTVIUI 550
#define VI_INTVIVI 551
#define VI_INTVIWI 552
#define VI_INTWIUI 553
#define VI_INTWIVI 554
#define VI_INTWIWI 555
#define WI_INTUIUI 556
#define WI_INTUIVI 557
#define WI_INTUIWI 558
#define WI_INTVIUI 559
#define WI_INTVIVI 560
#define WI_INTVIWI 561
#define WI_INTWIUI 562
#define WI_INTWIVI 563
#define WI_INTWIWI 564
#define VI_INTP_INT 565
#define WI_INTP_INT 566
#define dUINTdxdUINTdx 567
#define dUINTdxdUINTdy 568
#define dUINTdxdUINTdz 569
#define dUINTdydUINTdy 570
#define dUINTdydUINTdz 571
#define dUINTdzdUINTdz 572
#define dUINTdxdVINTdx 573
#define dUINTdxdVINTdy 574
#define dUINTdxdVINTdz 575
#define dUINTdydVINTdy 576
#define dUINTdydVINTdz 577
#define dUINTdzdVINTdz 578
#define dUINTdxdWINTdx 579
#define dUINTdxdWINTdy 580
#define dUINTdxdWINTdz 581
#define dUINTdydWINTdy 582
#define dUINTdydWINTdz 583
#define dUINTdzdWINTdz 584
#define dVINTdxdUINTdx 585
#define dVINTdxdUINTdy 586
#define dVINTdxdUINTdz 587
#define dVINTdydUINTdy 588
#define dVINTdydUINTdz 589
#define dVINTdzdUINTdz 590
#define dVINTdxdVINTdx 591
#define dVINTdxdVINTdy 592
#define dVINTdxdVINTdz 593
#define dVINTdydVINTdy 594
#define dVINTdydVINTdz 595
#define dVINTdzdVINTdz 596
#define dVINTdxdWINTdx 597
#define dVINTdxdWINTdy 598
#define dVINTdxdWINTdz 599
#define dVINTdydWINTdy 600
#define dVINTdydWINTdz 601
#define dVINTdzdWINTdz 602
#define dWINTdxdUINTdx 603
#define dWINTdxdUINTdy 604
#define dWINTdxdUINTdz 605
#define dWINTdydUINTdy 606
#define dWINTdydUINTdz 607
#define dWINTdzdUINTdz 608
#define dWINTdxdVINTdx 609
#define dWINTdxdVINTdy 610
#define dWINTdxdVINTdz 611
#define dWINTdydVINTdy 612
#define dWINTdydVINTdz 613
#define dWINTdzdVINTdz 614
#define dWINTdxdWINTdx 615
#define dWINTdxdWINTdy 616
#define dWINTdxdWINTdz 617
#define dWINTdydWINTdy 618
#define dWINTdydWINTdz 619
#define dWINTdzdWINTdz 620
#define P_INTdUINTdx 621
#define P_INTdUINTdy 622
#define P_INTdUINTdz 623
#define P_INTdVINTdx 624
#define P_INTdVINTdy 625
#define P_INTdVINTdz 626
#define P_INTdWINTdx 627
#define P_INTdWINTdy 628
#define P_INTdWINTdz 629
#define DISSIP_INT 630
//#define DISSIP_VAP_INT 793+4
//#define TRUE_DISSIP_INT 793+5  //
//#define TRUE_DISSIP_VAP_INT 793+6  //
#define P_NOPERTURBE 631
#define U_NOPERTURBE 632
#define V_NOPERTURBE 633
#define W_NOPERTURBE 634
#define I_NP 635
#define UI_INTP_INT 636
//// precaution pour le post de pression, permet de reperer d eventuel perturbation du champ de pression
//// toute les correlation contenant la pression sur l indicatrice non perturbee dans la phase liquide
#define P_NP 637
#define PIU_NP 638
#define PIV_NP 639
#define PIW_NP 640
#define IPdUdx_NP 641
#define IPdVdx_NP 642
#define IPdWdx_NP 643
#define IPdUdy_NP 644
#define IPdVdy_NP 645
#define IPdWdy_NP 646
#define IPdUdz_NP 647
#define IPdVdz_NP 648
#define IPdWdz_NP 649
#define IdPdx_NP 650
#define IdPdy_NP 651
#define IdPdz_NP 652
/// autre precaution : stat sur les pressions normale mais avec une valeur seuil a ne pas depasser.
/// Permet d'exclure de la moyenne les sauts de pression locaux trop important qui perturbe les stats
#define PI_seuil 653
#define PIU_seuil 654
#define PIV_seuil 655
#define PIW_seuil 656
#define IPdUdx_seuil 657
#define IPdVdx_seuil 658
#define IPdWdx_seuil 659
#define IPdUdy_seuil 660
#define IPdVdy_seuil 661
#define IPdWdy_seuil 662
#define IPdUdz_seuil 663
#define IPdVdz_seuil 664
#define IPdWdz_seuil 665
#define IdPdx_seuil 666
#define IdPdy_seuil 667
#define IdPdz_seuil 668
/// idem pour le champ de pression np
#define P_NP_seuil 669
#define PIU_NP_seuil 670
#define PIV_NP_seuil 671
#define PIW_NP_seuil 672
#define IPdUdx_NP_seuil 673
#define IPdVdx_NP_seuil 674
#define IPdWdx_NP_seuil 675
#define IPdUdy_NP_seuil 676
#define IPdVdy_NP_seuil 677
#define IPdWdy_NP_seuil 678
#define IPdUdz_NP_seuil 679
#define IPdVdz_NP_seuil 680
#define IPdWdz_NP_seuil 681
#define IdPdx_NP_seuil 682
#define IdPdy_NP_seuil 683
#define IdPdz_NP_seuil 684
#define Ivrappelx 685
#define Ivrappely 686
#define Ivrappelz 687
#define IP_INT 688
// stat pour travail turbulent de la force de rappel
#define UIvrappelx 689
#define VIvrappelx 690
#define WIvrappelx 691
#define UIvrappely 692
#define VIvrappely 693
#define WIvrappely 694
#define UIvrappelz 695
#define VIvrappelz 696
#define WIvrappelz 697

// Crossed-terms for the computation of the viscous stress tensor
#define dUdyaiNx_MOY 698
#define dUdzaiNx_MOY 699
#define dVdxaiNy_MOY 700
#define dVdzaiNy_MOY 701
#define dWdxaiNz_MOY 702
#define dWdyaiNz_MOY 703
//New field for the viscous stress in the vapor phase
#define dUdxIv_MOY 704
#define dUdyIv_MOY 705
#define dUdzIv_MOY 706
#define dVdxIv_MOY 707
#define dVdyIv_MOY 708
#define dVdzIv_MOY 709
#define dWdxIv_MOY 710
#define dWdyIv_MOY 711
#define dWdzIv_MOY 712
// Attempt to recompute the molecular diffusion in Rij
#define IddUdxdxU_MOY 713
#define IddVdxdxU_MOY 714
#define IddWdxdxU_MOY 715
#define IddUdxdxV_MOY 716
#define IddVdxdxV_MOY 717
#define IddWdxdxV_MOY 718
#define IddUdxdxW_MOY 719
#define IddVdxdxW_MOY 720
#define IddWdxdxW_MOY 721
#define IddUdydyU_MOY 722
#define IddVdydyU_MOY 723
#define IddWdydyU_MOY 724
#define IddUdydyV_MOY 725
#define IddVdydyV_MOY 726
#define IddWdydyV_MOY 727
#define IddUdydyW_MOY 728
#define IddVdydyW_MOY 729
#define IddWdydyW_MOY 730
#define IddUdzdzU_MOY 731
#define IddVdzdzU_MOY 732
#define IddWdzdzU_MOY 733
#define IddUdzdzV_MOY 734
#define IddVdzdzV_MOY 735
#define IddWdzdzV_MOY 736
#define IddUdzdzW_MOY 737
#define IddVdzdzW_MOY 738
#define IddWdzdzW_MOY 739
// Terms to reconstruct the physical pressure in the Rij equation In the post-processing they will need to be properly multiplyed by gravity constant and averaged density
// General average
#define Ix_MOY 740
#define xaiNx_MOY 741
#define xaiNy_MOY 742
#define xaiNz_MOY 743
// Pressure diffusion corrective term
#define UIx_MOY 744
#define VIx_MOY 745
#define WIx_MOY 746
// Redistribution corrective term
#define dUdxIx_MOY 747
#define dUdyIx_MOY 748
#define dUdzIx_MOY 749
#define dVdxIx_MOY 750
#define dVdyIx_MOY 751
#define dVdzIx_MOY 752
#define dWdxIx_MOY 753
#define dWdyIx_MOY 754
#define dWdzIx_MOY 755
// Interfacial terms
#define xUaiNx_MOY 756
#define xVaiNx_MOY 757
#define xWaiNx_MOY 758
#define xUaiNy_MOY 759
#define xVaiNy_MOY 760
#define xWaiNy_MOY 761
#define xUaiNz_MOY 762
#define xVaiNz_MOY 763
#define xWaiNz_MOY 764
// stat pour le champs de pression etendues
//ATT!!! nonostante si chiamino MOY, vi e comunque un if sui valori di pressione da considerare
#define BEG_VAP 765
#define P_VAP_Iv_MOY 765 //740
#define P_VAP_aiNx_MOY 766 //741
#define P_VAP_aiNy_MOY 767 //742
#define P_VAP_aiNz_MOY 768 //743
#define END_VAP 768
#define BEG_LIQ 769
// Momentum equation
#define P_LIQ_I_MOY 769 //744
#define P_LIQ_aiNx_MOY 770 //745
#define P_LIQ_aiNy_MOY 771 //746
#define P_LIQ_aiNz_MOY 772 //747
//Reynolds stress equation (only for liquid phase???)
#define UP_LIQ_aiNx_MOY 773 //748
#define VP_LIQ_aiNx_MOY 774 //749
#define WP_LIQ_aiNx_MOY 775 //750
#define UP_LIQ_aiNy_MOY 776 //751
#define VP_LIQ_aiNy_MOY 777 // 752
#define WP_LIQ_aiNy_MOY 778 //753
#define UP_LIQ_aiNz_MOY 779 //754
#define VP_LIQ_aiNz_MOY 780 //755
#define WP_LIQ_aiNz_MOY 781 //756
//Redistibution
#define IP_LIQ_dUdx_MOY 782 //757
#define IP_LIQ_dUdy_MOY 783 //758
#define IP_LIQ_dUdz_MOY 784 //759
#define IP_LIQ_dVdx_MOY 785 //760
#define IP_LIQ_dVdy_MOY 786 //761
#define IP_LIQ_dVdz_MOY 787 //762
#define IP_LIQ_dWdx_MOY 788 //763
#define IP_LIQ_dWdy_MOY 789 //764
#define IP_LIQ_dWdz_MOY 790 //765
//Diffusion
#define UP_LIQ_I_MOY 791 //766
#define VP_LIQ_I_MOY 792 //767
#define WP_LIQ_I_MOY 793 //768
#define END_LIQ 793

//#define BEG_DISSIP 793
#define DISSIP_VAP_MOY 794  // Après dernière ligne
#define TRUE_DISSIP_MOY 795
#define TRUE_DISSIP_VAP_MOY 796  // Dernière ligne
//#define END_DISSIP 796

//
