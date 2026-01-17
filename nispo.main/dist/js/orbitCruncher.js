(() => { // webpackBootstrap
"use strict";
var __webpack_modules__ = ({
"./src/engine/ootk/src/body/Celestial.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/body/Celestial.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Celestial: () => (Celestial)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Celestial is a static class that provides methods for calculating the position of celestial objects such as the Sun,
 * Moon, and planets in the sky. To create an instance of a Celestial object, use the Star class.
 */
class Celestial {
    constructor() {
        // disable constructor
    }
    /**
     * Calculates the azimuth and elevation of a celestial object at a given date, latitude,
     * longitude, right ascension, and declination.
     * @param date - The date for which to calculate the azimuth and elevation.
     * @param lat - The latitude of the observer.
     * @param lon - The longitude of the observer.
     * @param ra - The right ascension of the celestial object.
     * @param dec - The declination of the celestial object.
     * @returns An object containing the azimuth and elevation in degrees.
     */
    static azEl(date, lat, lon, ra, dec) {
        const c = {
            ra,
            dec,
            dist: 0,
        };
        const azEl = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.azEl(date, lat, lon, c);
        const el = (azEl.el + Celestial.atmosphericRefraction(azEl.el)); // elevation correction for refraction
        return {
            az: (azEl.az * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
            el: (el * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
        };
    }
    /**
     * Atmospheric refraction in astronomy, refers to the bending of light as it passes through the Earth's
     * atmosphere. This effect is most noticeable for celestial objects like stars and planets when they are
     * close to the horizon. Here's a breakdown of how it works:
     *
     * Actual Position: Due to this bending of light, the apparent position of a celestial object is slightly
     * different from its true position in the sky. When a star or planet is near the horizon, the effect is more
     * pronounced because the light path passes through more of the Earth's atmosphere, which increases the amount of
     * bending.
     *
     * A familiar example of atmospheric refraction is observed during sunrise and sunset. The Sun appears to
     * be above the horizon when it is actually just below it. This is because the light from the Sun is bent
     * upwards as it passes through the atmosphere.
     * @param h - elevation
     * @returns refraction
     */
    static atmosphericRefraction(h) {
        if (h < 0) {
            h = 0;
        }
        return (0.0002967 / Math.tan(h + 0.00312536 / (h + 0.08901179)));
    }
    /**
     * Calculate the declination. Similar to latitude on Earth, declination is another celestial coordinate.
     * It measures how far north or south an object is from the celestial equator
     * @param l - ecliptic longitude
     * @param b - ecliptic latitude
     * @returns declination
     */
    static declination(l, b) {
        return Math.asin(Math.sin(b) * Math.cos(_main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.e) + Math.cos(b) * Math.sin(_main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.e) * Math.sin(l));
    }
    /**
     * Calculate the right ascension. This is a celestial coordinate used to determine the position of objects
     * in the sky. It's analogous to longitude on Earth. Right Ascension indicates how far east an object is
     * from the vernal equinox along the celestial equator.
     * @param l - ecliptic longitude
     * @param b - ecliptic latitude
     * @returns right ascension
     */
    static rightAscension(l, b) {
        return Math.atan2(Math.sin(l) * Math.cos(_main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.e) - Math.tan(b) * Math.sin(_main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.e), Math.cos(l));
    }
    /**
     * Calculate the elevation. Elevation, or altitude, is the angle between an object in the sky and the
     * observer's local horizon. It's commonly expressed in degrees, where 0 degrees is right at the horizon
     * and 90 degrees is directly overhead (the zenith), but we are using radians to support trigonometric
     * functions like Math.sin() and Math.cos().
     * @param H - siderealTime
     * @param phi - latitude
     * @param dec - The declination of the sun
     * @returns elevation
     */
    static elevation(H, phi, dec) {
        return Math.asin(Math.sin(phi) * Math.sin(dec) + Math.cos(phi) * Math.cos(dec) * Math.cos(H));
    }
    /**
     * Calculate the azimuth. This is a compass direction measurement. Azimuth measures the angle along
     * the horizon from a specific reference direction (usually true north) to the point where a vertical
     * line from the object intersects the horizon.
     * @param H - siderealTime
     * @param phi - latitude
     * @param dec - The declination of the sun
     * @returns azimuth in rad
     */
    static azimuth(H, phi, dec) {
        return (Math.PI + Math.atan2(Math.sin(H), Math.cos(H) * Math.sin(phi) - Math.tan(dec) * Math.cos(phi)));
    }
}


}),
"./src/engine/ootk/src/body/Earth.ts": 
/*!*******************************************!*\
  !*** ./src/engine/ootk/src/body/Earth.ts ***!
  \*******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Earth: () => (Earth)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Earth metrics and operations.
class Earth {
    constructor() {
        // disable constructor
    }
    // / Earth gravitational parameter _(km²/s³)_.
    static mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.earthGravityParam;
    // / Earth equatorial radius.
    static radiusEquator = 6378.1363;
    // / Earth coefficient of flattening _(unitless)_.
    static flattening = 1.0 / 298.257223563;
    // / Earth polar radius.
    static radiusPolar = (Earth.radiusEquator * (1.0 - Earth.flattening));
    // / Earth mean radius.
    static radiusMean = ((2.0 * Earth.radiusEquator + Earth.radiusPolar) / 3.0);
    // / Earth eccentricity squared _(unitless)_.
    static eccentricitySquared = Earth.flattening * (2.0 - Earth.flattening);
    // / Earth J2 effect coefficient _(unitless)_.
    static j2 = 1.08262668355315e-3;
    // / Earth J3 effect coefficient _(unitless)_.
    static j3 = -2.53265648533224e-6;
    // / Earth J4 effect coefficient _(unitless)_.
    static j4 = -1.619621591367e-6;
    // / Earth J5 effect coefficient _(unitless)_.
    static j5 = -2.27296082868698e-7;
    // / Earth J6 effect coefficient _(unitless)_.
    static j6 = 5.40681239107085e-7;
    // / Earth rotation vector _(rad/s)_.
    static rotation = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(0, 0, 7.292115146706979e-5);
    // / Calculate mean motion _(rad/s)_ from a given [semimajorAxis] _(km)_.
    static smaToMeanMotion(semimajorAxis) {
        return Math.sqrt(Earth.mu / (semimajorAxis * semimajorAxis * semimajorAxis));
    }
    /**
     * Converts revolutions per day to semi-major axis.
     * @param rpd - The number of revolutions per day.
     * @returns The semi-major axis value.
     */
    static revsPerDayToSma(rpd) {
        return Earth.mu ** (1 / 3) / ((_main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * rpd) / _main_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerDay) ** (2 / 3);
    }
    // / Calculate Earth [PrecessionAngles] at a given UTC [epoch].
    static precession(epoch) {
        const t = epoch.toTT().toJulianCenturies();
        const zeta = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.zetaPoly_);
        const theta = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.thetaPoly_);
        const zed = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.zedPoly_);
        return { zeta: zeta, theta: theta, zed: zed };
    }
    // / Calculate Earth [NutationAngles] for a given UTC [epoch].
    static nutation(epoch) {
        const t = epoch.toTT().toJulianCenturies();
        const moonAnom = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.moonAnomPoly_);
        const sunAnom = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.sunAnomPoly_);
        const moonLat = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.moonLatPoly_);
        const sunElong = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.sunElongPoly_);
        const moonRaan = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.moonRaanPoly_);
        let deltaPsi = 0.0;
        let deltaEpsilon = 0.0;
        const dh = _main_js__WEBPACK_IMPORTED_MODULE_0__.DataHandler.getInstance();
        for (let i = 0; i < 4; i++) {
            const [a1, a2, a3, a4, a5, ai, bi, ci, di] = dh.getIau1980Coeffs(i);
            const arg = a1 * moonAnom + a2 * sunAnom + a3 * moonLat + a4 * sunElong + a5 * moonRaan;
            const sinC = ai + bi * t;
            const cosC = ci + di * t;
            deltaPsi += sinC * Math.sin(arg);
            deltaEpsilon += cosC * Math.cos(arg);
        }
        deltaPsi *= _main_js__WEBPACK_IMPORTED_MODULE_0__.ttasec2rad;
        deltaEpsilon *= _main_js__WEBPACK_IMPORTED_MODULE_0__.ttasec2rad;
        const meanEpsilon = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly)(t, Earth.meanEpsilonPoly_);
        const epsilon = meanEpsilon + deltaEpsilon;
        const eqEq = deltaPsi * Math.cos(meanEpsilon) +
            0.00264 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad * Math.sin(moonRaan) +
            0.000063 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad * Math.sin(2.0 * moonRaan);
        const gast = epoch.gmstAngle() + eqEq;
        return {
            dPsi: deltaPsi,
            dEps: deltaEpsilon,
            mEps: meanEpsilon,
            eps: epsilon,
            eqEq: eqEq,
            gast: gast,
        };
    }
    // / Convert a [semimajorAxis] _(km)_ to an eastward drift rate _(rad/day)_.
    static smaToDrift(semimajorAxis) {
        const t = (_main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * Math.sqrt(semimajorAxis ** 3 / Earth.mu)) / _main_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerSiderealDay;
        return (1.0 - t) * _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
    }
    // / Convert a [semimajorAxis] _(km)_ to an eastward drift rate _(°/day)_.
    static smaToDriftDegrees(semimajorAxis) {
        return Earth.smaToDrift(semimajorAxis) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    // / Convert an eastward [driftRate] _(rad/day)_ to a semimajor-axis _(km)_.
    static driftToSemimajorAxis(driftRate) {
        const t = (-driftRate / _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU + 1) * _main_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerSiderealDay;
        return ((Earth.mu * t * t) / (4 * Math.PI * Math.PI)) ** (1 / 3);
    }
    // / Convert an eastward [driftRate] _(°/day)_ to a semimajor-axis _(km)_.
    static driftDegreesToSma(driftRate) {
        return Earth.driftToSemimajorAxis(_main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * driftRate);
    }
    /**
     * Calculates the diameter of the Earth based on the satellite position.
     * @param satPos The position of the satellite.
     * @returns The diameter of the Earth.
     */
    static diameter(satPos) {
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.angularDiameter)(Earth.radiusEquator * 2, satPos.magnitude(), _main_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Sphere);
    }
    // / Earth precession `zeta` polynomial coefficients.
    static zetaPoly_ = Float64Array.from([
        0.017998 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        0.30188 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        2306.2181 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        0.0,
    ]);
    // / Earth precession `theta` polynomial coefficients.
    static thetaPoly_ = Float64Array.from([
        -0.041833 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        -0.42665 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        2004.3109 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        0.0,
    ]);
    // / Earth precession `zed` polynomial coefficients.
    static zedPoly_ = Float64Array.from([
        0.018203 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        1.09468 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        2306.2181 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        0,
    ]);
    static moonAnomPoly_ = Float64Array.from([
        1.4343e-5 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        0.0088553 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        (1325.0 * 360.0 + 198.8675605) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        134.96340251 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
    ]);
    // / Earth nutation Sun anomaly polynomial coefficients.
    static sunAnomPoly_ = Float64Array.from([
        3.8e-8 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        -0.0001537 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        (99.0 * 360.0 + 359.0502911) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        357.52910918 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
    ]);
    // / Earth nutation Moon latitude polynomial coefficients.
    static moonLatPoly_ = Float64Array.from([
        -2.88e-7 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        -0.003542 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        (1342.0 * 360.0 + 82.0174577) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        93.27209062 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
    ]);
    // / Earth nutation Sun elongation polynomial coefficients.
    static sunElongPoly_ = Float64Array.from([
        1.831e-6 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        -0.0017696 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        (1236.0 * 360.0 + 307.1114469) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        297.85019547 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
    ]);
    // / Earth nutation Moon right-ascension polynomial coefficients.
    static moonRaanPoly_ = Float64Array.from([
        2.139e-6 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        0.0020756 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        -(5.0 * 360.0 + 134.1361851) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
        125.04455501 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD,
    ]);
    // / Earth nutation mean epsilon polynomial coefficients.
    static meanEpsilonPoly_ = Float64Array.from([
        0.001813 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        -0.00059 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        -46.815 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
        84381.448 * _main_js__WEBPACK_IMPORTED_MODULE_0__.asec2rad,
    ]);
}


}),
"./src/engine/ootk/src/body/Moon.ts": 
/*!******************************************!*\
  !*** ./src/engine/ootk/src/body/Moon.ts ***!
  \******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Moon: () => (Moon)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../time/EpochUTC.js */ "./src/engine/ootk/src/time/EpochUTC.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _Earth_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _Sun_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./Sun.js */ "./src/engine/ootk/src/body/Sun.ts");
/**
 * @author Theodore Kruczek.
 * @license MIT
 * @copyright (c) 2022-2025 Theodore Kruczek Permission is
 * hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the
 * Software without restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @copyright (c) 2011-2015, Vladimir Agafonkin
 * @copyright (c) 2022 Robert Gester https://github.com/hypnos3/suncalc3
 * @see suncalc.LICENSE.md
 * Some of the math in this file was originally created by Vladimir Agafonkin.
 * Robert Gester's update was referenced for documentation. There were a couple
 * of bugs in both versions so there will be some differences if you are
 * migrating from either to this library.
 *
 * suncalc is a JavaScript library for calculating sun/moon position and light
 * phases. https://github.com/mourner/suncalc
 * It was reworked and enhanced by Robert Gester.
 *
 * The original suncalc is released under the terms of the BSD 2-Clause License.
 * @see http://aa.quae.nl/en/reken/hemelpositie.html
 * moon calculations are based on formulas from this website
 */







// / Moon metrics and operations.
class Moon {
    constructor() {
        // disable constructor
    }
    // / Moon gravitational parameter _(km³/s²)_.
    static mu = 4902.799;
    // / Moon equatorial radius _(km)_.
    static radiusEquator = 1738.0;
    // / Calculate the Moon's ECI position _(km)_ for a given UTC [epoch].
    static eci(epoch = _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_2__.EpochUTC.fromDateTime(new Date())) {
        const jc = epoch.toJulianCenturies();
        const dtr = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD;
        const lamEcl = 218.32 +
            481267.8813 * jc +
            6.29 * Math.sin((134.9 + 477198.85 * jc) * dtr) -
            1.27 * Math.sin((259.2 - 413335.38 * jc) * dtr) +
            0.66 * Math.sin((235.7 + 890534.23 * jc) * dtr) +
            0.21 * Math.sin((269.9 + 954397.7 * jc) * dtr) -
            0.19 * Math.sin((357.5 + 35999.05 * jc) * dtr) -
            0.11 * Math.sin((186.6 + 966404.05 * jc) * dtr);
        const phiEcl = 5.13 * Math.sin((93.3 + 483202.03 * jc) * dtr) +
            0.28 * Math.sin((228.2 + 960400.87 * jc) * dtr) -
            0.28 * Math.sin((318.3 + 6003.18 * jc) * dtr) -
            0.17 * Math.sin((217.6 - 407332.2 * jc) * dtr);
        const pllx = 0.9508 +
            0.0518 * Math.cos((134.9 + 477198.85 * jc) * dtr) +
            0.0095 * Math.cos((259.2 - 413335.38 * jc) * dtr) +
            0.0078 * Math.cos((235.7 + 890534.23 * jc) * dtr) +
            0.0028 * Math.cos((269.9 + 954397.7 * jc) * dtr);
        const obq = 23.439291 - 0.0130042 * jc;
        const rMag = 1 / Math.sin(pllx * dtr);
        const r = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(rMag * Math.cos(phiEcl * dtr) * Math.cos(lamEcl * dtr), rMag *
            (Math.cos(obq * dtr) * Math.cos(phiEcl * dtr) * Math.sin(lamEcl * dtr) -
                Math.sin(obq * dtr) * Math.sin(phiEcl * dtr)), rMag *
            (Math.sin(obq * dtr) * Math.cos(phiEcl * dtr) * Math.sin(lamEcl * dtr) +
                Math.cos(obq * dtr) * Math.sin(phiEcl * dtr)));
        const rMOD = r.scale(_Earth_js__WEBPACK_IMPORTED_MODULE_5__.Earth.radiusEquator);
        const p = _Earth_js__WEBPACK_IMPORTED_MODULE_5__.Earth.precession(epoch);
        return rMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
    }
    /**
     * Calculates the illumination of the Moon at a given epoch.
     * @param epoch - The epoch in UTC.
     * @param origin - The origin vector. Defaults to the origin vector if not provided.
     * @returns The illumination of the Moon, ranging from 0 to 1.
     */
    static illumination(epoch, origin) {
        const orig = origin ?? _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D.origin;
        const sunPos = _Sun_js__WEBPACK_IMPORTED_MODULE_6__.Sun.position(epoch).subtract(orig);
        const moonPos = this.eci(epoch).subtract(orig);
        const phaseAngle = sunPos.angle(moonPos);
        return 0.5 * (1 - Math.cos(phaseAngle));
    }
    /**
     * Calculates the diameter of the Moon.
     * @param obsPos - The position of the observer.
     * @param moonPos - The position of the Moon.
     * @returns The diameter of the Moon.
     */
    static diameter(obsPos, moonPos) {
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.angularDiameter)(this.radiusEquator * 2, obsPos.subtract(moonPos).magnitude(), _main_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Sphere);
    }
    /**
     * calculations for illumination parameters of the moon, based on
     * http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and Chapter 48 of "Astronomical Algorithms" 2nd
     * edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
     * @param date Date object or timestamp for calculating moon-illumination
     * @returns result object of moon-illumination
     */
    // eslint-disable-next-line max-statements
    static getMoonIllumination(date) {
        const dateValue = date instanceof Date ? date.getTime() : date;
        const lunarDaysMs = 2551442778; // The duration in days of a lunar cycle is 29.53058770576 days.
        const firstNewMoon2000 = 947178840000; // first newMoon in the year 2000 2000-01-06 18:14
        const dateObj = new Date(dateValue);
        const d = _Sun_js__WEBPACK_IMPORTED_MODULE_6__.Sun.date2jSince2000(dateObj);
        const s = _Sun_js__WEBPACK_IMPORTED_MODULE_6__.Sun.raDec(dateObj);
        const m = Moon.moonCoords(d);
        const sdist = 149598000; // distance from Earth to Sun in km
        const phi = Math.acos(Math.sin(s.dec) * Math.sin(m.dec) + Math.cos(s.dec) * Math.cos(m.dec) * Math.cos(s.ra - m.ra));
        const inc = Math.atan2(sdist * Math.sin(phi), m.dist - sdist * Math.cos(phi));
        const angle = Math.atan2(Math.cos(s.dec) * Math.sin(s.ra - m.ra), Math.sin(s.dec) * Math.cos(m.dec) - Math.cos(s.dec) * Math.sin(m.dec) * Math.cos(s.ra - m.ra));
        const phaseValue = 0.5 + (0.5 * inc * (angle < 0 ? -1 : 1)) / Math.PI;
        /*
         * calculates the difference in ms between the sirst fullMoon 2000 and given
         * Date
         */
        const diffBase = dateValue - firstNewMoon2000;
        // Calculate modulus to drop completed cycles
        let cycleModMs = diffBase % lunarDaysMs;
        // If negative number (date before new moon 2000) add lunarDaysMs
        if (cycleModMs < 0) {
            cycleModMs += lunarDaysMs;
        }
        const nextNewMoon = lunarDaysMs - cycleModMs + dateValue;
        let nextFullMoon = lunarDaysMs / 2 - cycleModMs + dateValue;
        if (nextFullMoon < dateValue) {
            nextFullMoon += lunarDaysMs;
        }
        const quater = lunarDaysMs / 4;
        let nextFirstQuarter = quater - cycleModMs + dateValue;
        if (nextFirstQuarter < dateValue) {
            nextFirstQuarter += lunarDaysMs;
        }
        let nextThirdQuarter = lunarDaysMs - quater - cycleModMs + dateValue;
        if (nextThirdQuarter < dateValue) {
            nextThirdQuarter += lunarDaysMs;
        }
        /*
         * Calculate the fraction of the moon cycle const currentfrac = cycleModMs /
         * lunarDaysMs;
         */
        const next = Math.min(nextNewMoon, nextFirstQuarter, nextFullMoon, nextThirdQuarter);
        // eslint-disable-next-line init-declarations
        let phase = null;
        for (const moonCycle of Moon.moonCycles_) {
            if (phaseValue >= moonCycle.from && phaseValue <= moonCycle.to) {
                phase = moonCycle;
                break;
            }
        }
        if (!phase) {
            throw new Error('Moon phase not found');
        }
        let type = '';
        if (next === nextNewMoon) {
            type = 'newMoon';
        }
        else if (next === nextFirstQuarter) {
            type = 'firstQuarter';
        }
        else if (next === nextFullMoon) {
            type = 'fullMoon';
        }
        else {
            type = 'thirdQuarter';
        }
        return {
            fraction: (1 + Math.cos(inc)) / 2,
            phase,
            phaseValue,
            angle,
            next: {
                value: next,
                date: new Date(next).toISOString(),
                type,
                newMoon: {
                    value: nextNewMoon,
                    date: new Date(nextNewMoon).toISOString(),
                },
                fullMoon: {
                    value: nextFullMoon,
                    date: new Date(nextFullMoon).toISOString(),
                },
                firstQuarter: {
                    value: nextFirstQuarter,
                    date: new Date(nextFirstQuarter).toISOString(),
                },
                thirdQuarter: {
                    value: nextThirdQuarter,
                    date: new Date(nextThirdQuarter).toISOString(),
                },
            },
        };
    }
    static rae(date, lat, lon) {
        const lw = (_utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * -lon);
        const phi = (_utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * lat);
        const d = _Sun_js__WEBPACK_IMPORTED_MODULE_6__.Sun.date2jSince2000(date);
        const c = Moon.moonCoords(d);
        const H = _Sun_js__WEBPACK_IMPORTED_MODULE_6__.Sun.siderealTime(d, lw) - c.ra;
        let h = _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.elevation(H, phi, c.dec);
        /*
         * formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus
         * (Willmann-Bell, Richmond) 1998.
         */
        const pa = Math.atan2(Math.sin(H), Math.tan(phi) * Math.cos(c.dec) - Math.sin(c.dec) * Math.cos(H));
        h = (h + _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.atmosphericRefraction(h)); // altitude correction for refraction
        return {
            az: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.azimuth(H, phi, c.dec),
            el: h,
            rng: c.dist,
            parallacticAngle: pa,
        };
    }
    /**
     * calculations for moon rise/set times are based on http://www.stargazing.net/kepler/moonrise.html article
     * @param date Date object or timestamp for calculating moon rise/set
     * @param lat Latitude of observer in degrees
     * @param lon Longitude of observer in degrees
     * @param isUtc If true, date will be interpreted as UTC
     * @returns result object of moon rise/set
     */
    static getMoonTimes(date, lat, lon, isUtc = false) {
        // Clone the date so we don't change the original
        const date_ = new Date(date);
        if (isUtc) {
            date_.setUTCHours(0, 0, 0, 0);
        }
        else {
            date_.setHours(0, 0, 0, 0);
        }
        const { rise, set, ye } = Moon.calculateRiseSetTimes_(date_, lat, lon);
        const result = {
            rise: null,
            set: null,
            ye: null,
            alwaysUp: null,
            alwaysDown: null,
            highest: null,
        };
        if (rise) {
            result.rise = new Date(Moon.hoursLater_(date_, rise));
        }
        if (set) {
            result.set = new Date(Moon.hoursLater_(date_, set));
        }
        if (!rise && !set) {
            if (ye > 0) {
                result.alwaysUp = true;
                result.alwaysDown = false;
            }
            else {
                result.alwaysUp = false;
                result.alwaysDown = true;
            }
        }
        else if (rise && set) {
            result.alwaysUp = false;
            result.alwaysDown = false;
            result.highest = new Date(Moon.hoursLater_(date_, Math.min(rise, set) + Math.abs(set - rise) / 2));
        }
        else {
            result.alwaysUp = false;
            result.alwaysDown = false;
        }
        return result;
    }
    static hoursLater_(date, h) {
        return new Date(date.getTime() + (h * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.MS_PER_DAY) / 24);
    }
    /**
     * Calculates the geocentric ecliptic coordinates of the moon.
     * @param d - The number of days since year 2000.
     * @returns An object containing the right ascension, declination, and
     * distance to the moon.
     */
    static moonCoords(d) {
        const L = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * (218.316 + 13.176396 * d); // ecliptic longitude
        const M = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * (134.963 + 13.064993 * d); // mean anomaly
        const F = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * (93.272 + 13.22935 * d); // mean distance
        const l = L + _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * 6.289 * Math.sin(M); // longitude
        const b = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD * 5.128 * Math.sin(F); // latitude
        const dt = 385001 - 20905 * Math.cos(M); // distance to the moon in km
        return {
            ra: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.rightAscension(l, b),
            dec: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.declination(l, b),
            dist: dt,
        };
    }
    static calculateRiseSetTimes_(t, lat, lon) {
        const hc = 0.133 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD;
        let h0 = Moon.rae(t, lat, lon).el - hc;
        let h1 = 0;
        let h2 = 0;
        let rise = 0;
        let set = 0;
        let a = 0;
        let b = 0;
        let xe = 0;
        let ye = 0;
        let d = 0;
        let roots = 0;
        let x1 = 0;
        let x2 = 0;
        let dx = 0;
        /*
         * go in 2-hour chunks, each time seeing if a 3-point quadratic curve
         * crosses zero (which means rise or set)
         */
        for (let i = 1; i <= 24; i += 2) {
            h1 = Moon.rae(Moon.hoursLater_(t, i), lat, lon).el - hc;
            h2 = Moon.rae(Moon.hoursLater_(t, i + 1), lat, lon).el - hc;
            a = (h0 + h2) / 2 - h1;
            b = (h2 - h0) / 2;
            xe = -b / (2 * a);
            ye = (a * xe + b) * xe + h1;
            d = b * b - 4 * a * h1;
            roots = 0;
            if (d >= 0) {
                dx = Math.sqrt(d) / (Math.abs(a) * 2);
                x1 = xe - dx;
                x2 = xe + dx;
                if (Math.abs(x1) <= 1) {
                    roots++;
                }
                if (Math.abs(x2) <= 1) {
                    roots++;
                }
                if (x1 < -1) {
                    x1 = x2;
                }
            }
            if (roots === 1) {
                if (h0 < 0) {
                    rise = i + x1;
                }
                else {
                    set = i + x1;
                }
            }
            else if (roots === 2) {
                rise = i + (ye < 0 ? x2 : x1);
                set = i + (ye < 0 ? x1 : x2);
            }
            if (rise && set) {
                break;
            }
            h0 = h2;
        }
        return { rise, set, ye };
    }
    static moonCycles_ = [
        {
            from: 0,
            to: 0.033863193308711,
            id: 'newMoon',
            emoji: '🌚',
            code: ':new_moon_with_face:',
            name: 'New Moon',
            weight: 1,
            css: 'wi-moon-new',
        },
        {
            from: 0.033863193308711,
            to: 0.216136806691289,
            id: 'waxingCrescentMoon',
            emoji: '🌒',
            code: ':waxing_crescent_moon:',
            name: 'Waxing Crescent',
            weight: 6.3825,
            css: 'wi-moon-wax-cres',
        },
        {
            from: 0.216136806691289,
            to: 0.283863193308711,
            id: 'firstQuarterMoon',
            emoji: '🌓',
            code: ':first_quarter_moon:',
            name: 'First Quarter',
            weight: 1,
            css: 'wi-moon-first-quart',
        },
        {
            from: 0.283863193308711,
            to: 0.466136806691289,
            id: 'waxingGibbousMoon',
            emoji: '🌔',
            code: ':waxing_gibbous_moon:',
            name: 'Waxing Gibbous',
            weight: 6.3825,
            css: 'wi-moon-wax-gibb',
        },
        {
            from: 0.466136806691289,
            to: 0.533863193308711,
            id: 'fullMoon',
            emoji: '🌝',
            code: ':full_moon_with_face:',
            name: 'Full Moon',
            weight: 1,
            css: 'wi-moon-full',
        },
        {
            from: 0.533863193308711,
            to: 0.716136806691289,
            id: 'waningGibbousMoon',
            emoji: '🌖',
            code: ':waning_gibbous_moon:',
            name: 'Waning Gibbous',
            weight: 6.3825,
            css: 'wi-moon-wan-gibb',
        },
        {
            from: 0.716136806691289,
            to: 0.783863193308711,
            id: 'thirdQuarterMoon',
            emoji: '🌗',
            code: ':last_quarter_moon:',
            name: 'third Quarter',
            weight: 1,
            css: 'wi-moon-third-quart',
        },
        {
            from: 0.783863193308711,
            to: 0.966136806691289,
            id: 'waningCrescentMoon',
            emoji: '🌘',
            code: ':waning_crescent_moon:',
            name: 'Waning Crescent',
            weight: 6.3825,
            css: 'wi-moon-wan-cres',
        },
        {
            from: 0.966136806691289,
            to: 1,
            id: 'newMoon',
            emoji: '🌚',
            code: ':new_moon_with_face:',
            name: 'New Moon',
            weight: 1,
            css: 'wi-moon-new',
        },
    ];
}


}),
"./src/engine/ootk/src/body/Sun.ts": 
/*!*****************************************!*\
  !*** ./src/engine/ootk/src/body/Sun.ts ***!
  \*****************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sun: () => (Sun)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 * @copyright (c) 2011-2015, Vladimir Agafonkin
 * @copyright (c) 2022 Robert Gester https://github.com/hypnos3/suncalc3
 * @see suncalc.LICENSE.md
 * Some of the math in this file was originally created by Vladimir Agafonkin.
 * Robert Gester's update was referenced for documentation. There were a couple
 * of bugs in both versions so there will be some differences if you are
 * migrating from either to this library.
 *
 * suncalc is a JavaScript library for calculating sun/moon position and light
 * phases. https://github.com/mourner/suncalc
 * It was reworked and enhanced by Robert Gester.
 *
 * The original suncalc is released under the terms of the BSD 2-Clause License.
 */

/**
 * Sun metrics and operations.
 */
class Sun {
    static J0_ = 0.0009;
    static J1970_ = 2440587.5;
    static J2000_ = 2451545;
    static e = _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * 23.4397;
    /**
     * Array representing the times for different phases of the sun. Each element
     * in the array contains:
     * - The angle in degrees representing the time offset from solar noon.
     * - The name of the start time for the phase.
     * - The name of the end time for the phase.
     */
    static times_ = [
        [6, 'goldenHourDawnEnd', 'goldenHourDuskStart'], // GOLDEN_HOUR_2
        [-0.3, 'sunriseEnd', 'sunsetStart'], // SUNRISE_END
        [-0.833, 'sunriseStart', 'sunsetEnd'], // SUNRISE
        [-1, 'goldenHourDawnStart', 'goldenHourDuskEnd'], // GOLDEN_HOUR_1
        [-4, 'blueHourDawnEnd', 'blueHourDuskStart'], // BLUE_HOUR
        [-6, 'civilDawn', 'civilDusk'], // DAWN
        [-8, 'blueHourDawnStart', 'blueHourDuskEnd'], // BLUE_HOUR
        [-12, 'nauticalDawn', 'nauticalDusk'], // NAUTIC_DAWN
        [-15, 'amateurDawn', 'amateurDusk'],
        [-18, 'astronomicalDawn', 'astronomicalDusk'], // ASTRO_DAWN
    ];
    /**
     * Gravitational parameter of the Sun. (km³/s²)
     */
    static mu = 1.32712428e11;
    /**
     * The angle of the penumbra of the Sun, in radians.
     */
    static penumbraAngle = (0.26900424 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    /**
     * The radius of the Sun in kilometers.
     */
    static radius = 695500.0;
    /**
     * The mean solar flux of the Sun. (W/m²)
     */
    static solarFlux = 1367.0;
    /**
     * The solar pressure exerted by the Sun. (N/m²) It is calculated as the solar
     * flux divided by the speed of light.
     */
    static solarPressure = Sun.solarFlux / (_main_js__WEBPACK_IMPORTED_MODULE_0__.cKmPerSec * 1000);
    /**
     * The angle of the umbra, in radians.
     */
    static umbraAngle = (0.26411888 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    constructor() {
        // disable constructor
    }
    /**
     * Calculates the azimuth and elevation of the Sun for a given date, latitude,
     * and longitude.
     * @param date - The date for which to calculate the azimuth and elevation.
     * @param lat - The latitude in degrees.
     * @param lon - The longitude in degrees.
     * @param c - The right ascension and declination of the target. Defaults to
     * the Sun's right ascension and declination
     * @returns An object containing the azimuth and elevation of the Sun in
     * radians.
     */
    static azEl(date, lat, lon, c) {
        const lw = (-lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const phi = (lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const d = Sun.date2jSince2000(date);
        c ??= Sun.raDec(date);
        const H = Sun.siderealTime(d, lw) - c.ra;
        return {
            az: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.azimuth(H, phi, c.dec),
            el: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.elevation(H, phi, c.dec),
        };
    }
    /**
     * get number of days for a dateValue since 2000
     * See: https://en.wikipedia.org/wiki/Epoch_(astronomy)
     * @param date date as timestamp to get days
     * @returns count of days
     */
    static date2jSince2000(date) {
        return date.getTime() / _main_js__WEBPACK_IMPORTED_MODULE_0__.MS_PER_DAY + Sun.J1970_ - Sun.J2000_;
    }
    /**
     * Calculates the angular diameter of the Sun given the observer's position
     * and the Sun's position.
     * @param obsPos The observer's position in kilometers.
     * @param sunPos The Sun's position in kilometers.
     * @returns The angular diameter of the Sun in radians.
     */
    static diameter(obsPos, sunPos) {
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.angularDiameter)(this.radius * 2, obsPos.subtract(sunPos).magnitude(), _main_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Sphere);
    }
    /**
     * Calculate eclipse angles given a satellite ECI position and Sun ECI
     * position.
     * @param satPos the satellite position
     * @param sunPos the sun position
     * @returns [central body angle, central body apparent radius, sun apparent]
     */
    static eclipseAngles(satPos, sunPos) {
        const satSun = sunPos.subtract(satPos);
        const r = satPos.magnitude();
        return [
            // central body angle
            satSun.angle(satPos.negate()),
            // central body apparent radius
            Math.asin(_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator / r),
            // sun apparent radius
            Math.asin(this.radius / satSun.magnitude()),
        ];
    }
    /**
     * Ecliptic latitude measures the distance north or south of the ecliptic,
     * attaining +90° at the north ecliptic pole (NEP) and -90° at the south
     * ecliptic pole (SEP). The ecliptic itself is 0° latitude.
     * @param B - ?
     * @returns ecliptic latitude
     */
    static eclipticLatitude(B) {
        const C = _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU / 360;
        const L = B - 0.00569 - 0.00478 * Math.sin(C * B);
        return _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * (L + 0.0003 * Math.sin(C * 2 * L));
    }
    /**
     * Ecliptic longitude, also known as celestial longitude, measures the angular
     * distance of an object along the ecliptic from the primary direction. It is
     * measured positive eastwards in the fundamental plane (the ecliptic) from 0°
     * to 360°. The primary direction (0° ecliptic longitude) points from the
     * Earth towards the Sun at the vernal equinox of the Northern Hemisphere. Due
     * to axial precession, the ecliptic longitude of most "fixed stars" increases
     * by about 50.3 arcseconds per year, or 83.8 arcminutes per century.
     * @param M - solar mean anomaly
     * @returns ecliptic longitude
     */
    static eclipticLongitude(M) {
        const C = _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * (1.9148 * Math.sin(M) + 0.02 * Math.sin(2 * M) + 0.0003 * Math.sin(3 * M));
        const P = _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * 102.9372; // perihelion of Earth
        return (M + C + P + Math.PI); // Sun's mean longitude
    }
    /**
     * returns set time for the given sun altitude
     * @param h - height at 0
     * @param lw - rad * -lng
     * @param phi -  rad * lat;
     * @param dec - declination
     * @param n - Julian cycle
     * @param M - solar mean anomal
     * @param L - ecliptic longitude
     * @returns set time
     */
    static getSetJulian(h, lw, phi, dec, n, M, L) {
        const w = Sun.hourAngle(h, phi, dec);
        const a = Sun.approxTransit_(w, lw, n);
        return Sun.solarTransitJulian(a, M, L);
    }
    /**
     * Calculates the time of the sun based on the given azimuth.
     * @param dateValue - The date value or Date object.
     * @param lat - The latitude in degrees.
     * @param lon - The longitude in degrees.
     * @param az - The azimuth in radians or degrees.
     * @param isDegrees - Indicates if the azimuth is in degrees. Default is false.
     * @returns The calculated time of the sun.
     * @throws Error if the azimuth, latitude, or longitude is missing.
     */
    static getSunTimeByAz(dateValue, lat, lon, az, isDegrees = false) {
        if (isNaN(az)) {
            throw new Error('azimuth missing');
        }
        if (isNaN(lat)) {
            throw new Error('latitude missing');
        }
        if (isNaN(lon)) {
            throw new Error('longitude missing');
        }
        if (isDegrees) {
            az = (az * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        }
        const date = dateValue instanceof Date ? dateValue : new Date(dateValue);
        const lw = (_main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * -lon);
        const phi = (_main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * lat);
        let dateVal = new Date(date.getFullYear(), date.getMonth(), date.getDate(), 0, 0, 0).getTime();
        let addval = _main_js__WEBPACK_IMPORTED_MODULE_0__.MS_PER_DAY; // / 2);
        dateVal += addval;
        while (addval > 200) {
            const newDate = new Date(dateVal);
            const d = Sun.date2jSince2000(newDate);
            const c = Sun.raDec(newDate);
            const H = Sun.siderealTime(d, lw) - c.ra;
            const newAz = _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.azimuth(H, phi, c.dec);
            addval /= 2;
            if (newAz < az) {
                dateVal += addval;
            }
            else {
                dateVal -= addval;
            }
        }
        return new Date(Math.floor(dateVal));
    }
    /**
     * Calculates sun times for a given date and latitude/longitude
     *
     * Default altitude is 0 meters. If `isUtc` is `true`, the times are returned
     * as UTC, otherwise in local time.
     * @param dateVal - The date value or Date object.
     * @param lat - The latitude in degrees.
     * @param lon - The longitude in degrees.
     * @param alt - The altitude in meters. Default is 0.
     * @param isUtc - Indicates if the times should be returned as UTC. Default is
     * false.
     * @returns An object containing the times of the sun.
     */
    static getTimes(dateVal, lat, lon, alt = 0, isUtc = false) {
        if (isNaN(lat)) {
            throw new Error('latitude missing');
        }
        if (isNaN(lon)) {
            throw new Error('longitude missing');
        }
        // Ensure date is a Date object
        const date = dateVal instanceof Date ? dateVal : new Date(dateVal);
        if (isUtc) {
            date.setUTCHours(12, 0, 0, 0);
        }
        else {
            date.setHours(12, 0, 0, 0);
        }
        let time;
        let h0 = 0;
        let Jset = 0;
        let Jrise = 0;
        const { Jnoon, dh, lw, phi, dec, n, M, L } = Sun.calculateJnoon_(lon, lat, alt, date);
        // Determine when the sun is at its highest and lowest (darkest) points.
        const result = {
            solarNoon: Sun.julian2date(Jnoon),
            nadir: Sun.julian2date(Jnoon + 0.5), // https://github.com/mourner/suncalc/pull/125
        };
        // Add all other unique times using Jnoon as a reference
        for (let i = 0, len = Sun.times_.length; i < len; i += 1) {
            time = Sun.times_[i];
            const angle = time[0];
            h0 = ((angle + dh) * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
            Jset = Sun.getSetJ_(h0, lw, phi, dec, n, M, L);
            Jrise = Jnoon - (Jset - Jnoon);
            result[time[1]] = Sun.julian2date(Jrise);
            result[time[2]] = Sun.julian2date(Jset);
        }
        return result;
    }
    /**
     * hour angle
     * @param h - heigh at 0
     * @param phi -  rad * lat;
     * @param dec - declination
     * @returns hour angle
     */
    static hourAngle(h, phi, dec) {
        return Math.acos((Math.sin(h) - Math.sin(phi) * Math.sin(dec)) / (Math.cos(phi) * Math.cos(dec)));
    }
    /**
     * convert Julian calendar to date object
     * @param julian day number in Julian calendar to convert
     * @returns result date as timestamp
     */
    static julian2date(julian) {
        return new Date((julian - Sun.J1970_) * _main_js__WEBPACK_IMPORTED_MODULE_0__.MS_PER_DAY);
    }
    /**
     * Julian cycle
     *
     * The Julian cycle is a period of 7980 years after which the positions of the
     * Sun, Moon, and planets repeat. It is used in astronomical calculations to
     * determine the position of celestial bodies.
     *
     * The Julian Period starts at noon on January 1, 4713 B.C.E. (Julian
     * calendar) and lasts for 7980 years. This was determined because it is a
     * time period long enough to include all of recorded history and includes
     * some time in the future that would incorporate the three important
     * calendrical cycles, the Golden Number Cycle, the Solar Cycle, and the Roman
     * Indiction.
     *
     * The Golden Number Cycle is a cycle of 19 years, while the Solar Cycle is a
     * cycle of 28 years and the Roman Indiction repeats every 15 years. Thus the
     * Julian Period is calculated to be 7980 years long or 2,914,695 days because
     * 19*28*15 = 7980.
     * @param date - Date object for calculating julian cycle
     * @param lon - Degrees longitude
     * @returns julian cycle
     */
    static julianCycle(date, lon) {
        const lw = (-lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const d = Sun.date2jSince2000(date);
        return Math.round(d - Sun.J0_ - lw / ((2 * _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU) / 2));
    }
    /**
     * Calculate the lighting ratio given a satellite ECI position [satPos] _(km)_
     * and Sun ECI position [sunPos] _(km)_.
     *
     * Returns `1.0` if the satellite is fully illuminated and `0.0` when fully
     * eclipsed.
     * @param satPos - The position of the satellite.
     * @param sunPos - The position of the sun.
     * @returns The lighting ratio.
     */
    static lightingRatio(satPos, sunPos) {
        const [sunSatAngle, aCent, aSun] = Sun.eclipseAngles(satPos, sunPos);
        if (sunSatAngle - aCent + aSun <= 1e-10) {
            return 0.0;
        }
        else if (sunSatAngle - aCent - aSun < -1e-10) {
            const ssa2 = sunSatAngle * sunSatAngle;
            const ssaInv = 1.0 / (2.0 * sunSatAngle);
            const ac2 = aCent * aCent;
            const as2 = aSun * aSun;
            const acAsDiff = ac2 - as2;
            const a1 = (ssa2 - acAsDiff) * ssaInv;
            const a2 = (ssa2 + acAsDiff) * ssaInv;
            const asr1 = a1 / aSun;
            const asr2 = as2 - a1 * a1;
            const acr1 = a2 / aCent;
            const acr2 = ac2 - a2 * a2;
            const p1 = as2 * Math.acos(asr1) - a1 * Math.sqrt(asr2);
            const p2 = ac2 * Math.acos(acr1) - a2 * Math.sqrt(acr2);
            return 1.0 - (p1 + p2) / (Math.PI * as2);
        }
        return 1.0;
    }
    /**
     * Calculates the lighting factor based on the position of the satellite and the sun.
     * @deprecated This method was previously used. It is now deprecated and will be removed
     * in a future release.
     * @param satPos The position of the satellite.
     * @param sunPos The position of the sun.
     * @returns The lighting factor.
     */
    static sunlightLegacy(satPos, sunPos) {
        let lighting = 1.0;
        const semiDiamEarth = Math.asin(_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean / Math.sqrt((-satPos.x) ** 2 + (-satPos.y) ** 2 + (-satPos.z) ** 2)) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
        const semiDiamSun = Math.asin(Sun.radius / Math.sqrt((-satPos.x + sunPos.x) ** 2 + (-satPos.y + sunPos.y) ** 2 + (-satPos.z + sunPos.z) ** 2)) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
        // Angle between earth and sun
        const theta = Math.acos(satPos.negate().dot(sunPos.negate()) /
            (Math.sqrt((-satPos.x) ** 2 + (-satPos.y) ** 2 + (-satPos.z) ** 2) *
                Math.sqrt((-satPos.x + sunPos.x) ** 2 + (-satPos.y + sunPos.y) ** 2 + (-satPos.z + sunPos.z) ** 2))) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
        if (semiDiamEarth > semiDiamSun && theta < semiDiamEarth - semiDiamSun) {
            lighting = 0;
        }
        if (Math.abs(semiDiamEarth - semiDiamSun) < theta && theta < semiDiamEarth + semiDiamSun) {
            lighting = 0.5;
        }
        if (semiDiamSun > semiDiamEarth) {
            lighting = 0.5;
        }
        if (theta < semiDiamSun - semiDiamEarth) {
            lighting = 0.5;
        }
        return lighting;
    }
    /**
     * Calculates the position vector of the Sun at a given epoch in the
     * Earth-centered inertial (ECI) coordinate system.
     * @param epoch - The epoch in UTC.
     * @returns The position vector of the Sun in Kilometers.
     */
    static position(epoch) {
        const jc = epoch.toJulianCenturies();
        const dtr = _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
        const lamSun = 280.46 + 36000.77 * jc;
        const mSun = 357.5291092 + 35999.05034 * jc;
        const lamEc = lamSun + 1.914666471 * Math.sin(mSun * dtr) + 0.019994643 * Math.sin(2.0 * mSun * dtr);
        const obliq = 23.439291 - 0.0130042 * jc;
        const rMag = 1.000140612 - 0.016708617 * Math.cos(mSun * dtr) - 0.000139589 * Math.cos(2.0 * mSun * dtr);
        const r = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(rMag * Math.cos(lamEc * dtr), rMag * Math.cos(obliq * dtr) * Math.sin(lamEc * dtr), rMag * Math.sin(obliq * dtr) * Math.sin(lamEc * dtr));
        const rMOD = r.scale(_main_js__WEBPACK_IMPORTED_MODULE_0__.astronomicalUnit);
        const p = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.precession(epoch);
        return rMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
    }
    /**
     * Calculate the Sun's apparent ECI position _(km)_ from Earth for a given UTC
     * [epoch].
     * @param epoch - The epoch in UTC.
     * @returns The Sun's apparent ECI position in kilometers.
     */
    static positionApparent(epoch) {
        const distance = Sun.position(epoch).magnitude();
        const dSec = distance / _main_js__WEBPACK_IMPORTED_MODULE_0__.cKmPerSec;
        return Sun.position(epoch.roll(-dSec));
    }
    /**
     * Calculates the right ascension and declination of the Sun for a given date.
     * @param date - The date for which to calculate the right ascension and declination.
     * @returns An object containing the declination and right ascension of the Sun.
     */
    static raDec(date) {
        const d = Sun.date2jSince2000(date);
        const M = Sun.solarMeanAnomaly_(d);
        const L = Sun.eclipticLongitude(M);
        return {
            dec: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.declination(L, 0),
            ra: _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.rightAscension(L, 0),
            dist: 0,
        };
    }
    /**
     * Return `true` if the ECI satellite position [posSat] is in eclipse at the
     * given UTC [epoch].
     * @param epoch - The epoch in UTC.
     * @param posSat - The ECI position of the satellite in kilometers.
     * @returns `true` if the satellite is in eclipse.
     */
    static shadow(epoch, posSat) {
        const posSun = Sun.positionApparent(epoch);
        let shadow = false;
        if (posSun.dot(posSat) < 0) {
            const angle = posSun.angle(posSat);
            const r = posSat.magnitude();
            const satHoriz = r * Math.cos(angle);
            const satVert = r * Math.sin(angle);
            const penVert = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator + Math.tan(this.penumbraAngle) * satHoriz;
            if (satVert <= penVert) {
                shadow = true;
            }
        }
        return shadow;
    }
    /**
     * side real time
     * @param d - julian day
     * @param lw - longitude of the observer
     * @returns sidereal time
     */
    static siderealTime(d, lw) {
        return _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * (280.16 + 360.9856235 * d) - lw;
    }
    /**
     * solar transit in Julian
     * @param ds approxTransit
     * @param M solar mean anomal
     * @param L ecliptic longitude
     * @returns solar transit in Julian
     */
    static solarTransitJulian(ds, M, L) {
        return Sun.J2000_ + ds + 0.0053 * Math.sin(M) - 0.0069 * Math.sin(2 * L);
    }
    /**
     * The approximate transit time
     * @param Ht hourAngle
     * @param lw rad * -lng
     * @param n Julian cycle
     * @returns approx transit
     */
    static approxTransit_(Ht, lw, n) {
        return Sun.J0_ + (Ht + lw) / _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU + n;
    }
    static calculateJnoon_(lon, lat, alt, date) {
        const lw = (_main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * -lon);
        const phi = (_main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * lat);
        const dh = Sun.observerAngle_(alt);
        const d = Sun.date2jSince2000(date);
        const n = Sun.julianCycle_(d, lw);
        const ds = Sun.approxTransit_(0, lw, n);
        const M = Sun.solarMeanAnomaly_(ds);
        const L = Sun.eclipticLongitude(M);
        const dec = _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.declination(L, 0);
        const Jnoon = Sun.solarTransitJulian(ds, M, L);
        return { Jnoon, dh, lw, phi, dec, n, M, L };
    }
    /**
     * returns set time for the given sun altitude
     * @param alt altitude at 0
     * @param lw lng
     * @param phi lat
     * @param dec declination
     * @param n Julian cycle
     * @param M solar mean anomal
     * @param L ecliptic longitude
     * @returns sunset time in days since 2000
     */
    static getSetJ_(alt, lw, phi, dec, n, M, L) {
        const w = Sun.hourAngle(alt, phi, dec);
        const a = Sun.approxTransit_(w, lw, n);
        return Sun.solarTransitJulian(a, M, L);
    }
    static julianCycle_(d, lw) {
        const lonOffset = lw / _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
        return Math.round(d - Sun.J0_ - lonOffset);
    }
    /**
     * calculates the obderver angle
     * @param alt the observer altitude (in meters) relative to the horizon
     * @returns height for further calculations
     */
    static observerAngle_(alt) {
        return ((-2.076 * Math.sqrt(alt)) / 60);
    }
    /**
     * get solar mean anomaly
     * @param d julian day
     * @returns solar mean anomaly
     */
    static solarMeanAnomaly_(d) {
        return _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD * (357.5291 + 0.98560028 * d);
    }
}


}),
"./src/engine/ootk/src/body/index.ts": 
/*!*******************************************!*\
  !*** ./src/engine/ootk/src/body/index.ts ***!
  \*******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Celestial: () => (/* reexport safe */ _Celestial_js__WEBPACK_IMPORTED_MODULE_0__.Celestial),
  Earth: () => (/* reexport safe */ _Earth_js__WEBPACK_IMPORTED_MODULE_1__.Earth),
  Moon: () => (/* reexport safe */ _Moon_js__WEBPACK_IMPORTED_MODULE_2__.Moon),
  Sun: () => (/* reexport safe */ _Sun_js__WEBPACK_IMPORTED_MODULE_3__.Sun)
});
/* ESM import */var _Celestial_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Celestial.js */ "./src/engine/ootk/src/body/Celestial.ts");
/* ESM import */var _Earth_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _Moon_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Moon.js */ "./src/engine/ootk/src/body/Moon.ts");
/* ESM import */var _Sun_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Sun.js */ "./src/engine/ootk/src/body/Sun.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






}),
"./src/engine/ootk/src/coordinate/ClassicalElements.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/ClassicalElements.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ClassicalElements: () => (ClassicalElements)
});
/* ESM import */var _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/OrbitRegime.js */ "./src/engine/ootk/src/enums/OrbitRegime.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _EquinoctialElements_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./EquinoctialElements.js */ "./src/engine/ootk/src/coordinate/EquinoctialElements.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






/**
 * The ClassicalElements class represents the classical orbital elements of an object.
 * @example
 * ```ts
 * const epoch = EpochUTC.fromDateTime(new Date('2024-01-14T14:39:39.914Z'));
 * const elements = new ClassicalElements({
 *  epoch,
 *  semimajorAxis: 6943.547853722985 as Kilometers,
 *  eccentricity: 0.0011235968124658146,
 *  inclination: 0.7509087232045765 as Radians,
 *  rightAscension: 0.028239555738616327 as Radians,
 *  argPerigee: 2.5386411901807353 as Radians,
 *  trueAnomaly: 0.5931399364974058 as Radians,
 * });
 * ```
 */
class ClassicalElements {
    epoch;
    semimajorAxis;
    eccentricity;
    inclination;
    rightAscension;
    argPerigee;
    trueAnomaly;
    /** Gravitational parameter in km³/s².  */
    mu;
    constructor({ epoch, semimajorAxis, eccentricity, inclination, rightAscension, argPerigee, trueAnomaly, mu = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.earthGravityParam, }) {
        this.epoch = epoch;
        this.semimajorAxis = semimajorAxis;
        this.eccentricity = eccentricity;
        this.inclination = inclination;
        this.rightAscension = rightAscension;
        this.argPerigee = argPerigee;
        this.trueAnomaly = trueAnomaly;
        this.mu = mu;
    }
    /**
     * Creates a new instance of ClassicalElements from a StateVector.
     * @param state The StateVector to convert.
     * @param mu The gravitational parameter of the central body. Default value is Earth's gravitational parameter.
     * @returns A new instance of ClassicalElements.
     * @throws Error if the StateVector is not in an inertial frame.
     */
    static fromStateVector(state, mu = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.earthGravityParam) {
        if (!state.inertial) {
            throw new Error('State vector must be in inertial frame (like J2000).');
        }
        const pos = state.position;
        const vel = state.velocity;
        const a = state.semimajorAxis;
        const eVecA = pos.scale(vel.magnitude() ** 2 - mu / pos.magnitude());
        const eVecB = vel.scale(pos.dot(vel));
        const eVec = eVecA.subtract(eVecB).scale(1 / mu);
        const e = eVec.magnitude();
        const h = pos.cross(vel);
        const i = Math.acos((0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.clamp)(h.z / h.magnitude(), -1.0, 1.0));
        const n = _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D.zAxis.cross(h);
        let o = Math.acos((0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.clamp)(n.x / n.magnitude(), -1.0, 1.0));
        if (n.y < 0) {
            o = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU - o;
        }
        let w = n.angle(eVec);
        if (eVec.z < 0) {
            w = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU - w;
        }
        let v = eVec.angle(pos);
        if (pos.dot(vel) < 0) {
            v = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU - v;
        }
        return new ClassicalElements({
            epoch: state.epoch,
            semimajorAxis: a,
            eccentricity: e,
            inclination: i,
            rightAscension: o,
            argPerigee: w,
            trueAnomaly: v,
            mu,
        });
    }
    /**
     * Gets the inclination in degrees.
     * @returns The inclination in degrees.
     */
    get inclinationDegrees() {
        return (this.inclination * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG);
    }
    /**
     * Gets the right ascension in degrees.
     * @returns The right ascension in degrees.
     */
    get rightAscensionDegrees() {
        return (this.rightAscension * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG);
    }
    /**
     * Gets the argument of perigee in degrees.
     * @returns The argument of perigee in degrees.
     */
    get argPerigeeDegrees() {
        return (this.argPerigee * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG);
    }
    /**
     * Gets the true anomaly in degrees.
     * @returns The true anomaly in degrees.
     */
    get trueAnomalyDegrees() {
        return (this.trueAnomaly * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG);
    }
    /**
     * Gets the apogee of the classical elements. It is measured from the surface of the earth.
     * @returns The apogee in kilometers.
     */
    get apogee() {
        return (this.semimajorAxis * (1.0 + this.eccentricity) - _main_js__WEBPACK_IMPORTED_MODULE_1__.Earth.radiusMean);
    }
    /**
     * Gets the perigee of the classical elements. The perigee is the point in an
     * orbit that is closest to the surface of the earth.
     * @returns The perigee distance in kilometers.
     */
    get perigee() {
        return (this.semimajorAxis * (1.0 - this.eccentricity) - _main_js__WEBPACK_IMPORTED_MODULE_1__.Earth.radiusMean);
    }
    toString() {
        return [
            '[ClassicalElements]',
            `  Epoch: ${this.epoch}`,
            `  Semimajor Axis (a):       ${this.semimajorAxis.toFixed(4)} km`,
            `  Eccentricity (e):         ${this.eccentricity.toFixed(7)}`,
            `  Inclination (i):          ${this.inclinationDegrees.toFixed(4)}°`,
            `  Right Ascension (Ω):      ${this.rightAscensionDegrees.toFixed(4)}°`,
            `  Argument of Perigee (ω):  ${this.argPerigeeDegrees.toFixed(4)}°`,
            `  True Anomaly (ν):         ${this.trueAnomalyDegrees.toFixed(4)}°`,
        ].join('\n');
    }
    /**
     * Calculates the mean motion of the celestial object.
     * @returns The mean motion in radians.
     */
    get meanMotion() {
        return Math.sqrt(this.mu / this.semimajorAxis ** 3);
    }
    /**
     * Calculates the period of the orbit.
     * @returns The period in seconds.
     */
    get period() {
        const periodSec = (_utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU * Math.sqrt(this.semimajorAxis ** 3 / this.mu));
        return (periodSec / 60);
    }
    /**
     * Compute the number of revolutions completed per day for this orbit.
     * @returns The number of revolutions per day.
     */
    get revsPerDay() {
        return _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.MINUTES_PER_DAY / this.period;
    }
    /**
     * Returns the orbit regime based on the classical elements.
     * @returns The orbit regime.
     */
    getOrbitRegime() {
        const n = this.revsPerDay;
        const p = this.period * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.sec2min;
        if (n >= 0.99 && n <= 1.01 && this.eccentricity < 0.01) {
            return _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime.GEO;
        }
        if (p >= 600 && p <= 800 && this.eccentricity <= 0.25) {
            return _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime.MEO;
        }
        if (n >= 11.25 && this.eccentricity <= 0.25) {
            return _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime.LEO;
        }
        if (this.eccentricity > 0.25) {
            return _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime.HEO;
        }
        return _enums_OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime.OTHER;
    }
    /**
     * Converts the classical orbital elements to position and velocity vectors.
     * @returns An object containing the position and velocity vectors.
     */
    toPositionVelocity() {
        const rVec = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D(Math.cos(this.trueAnomaly), Math.sin(this.trueAnomaly), 0.0);
        const rPQW = rVec.scale((this.semimajorAxis * (1.0 - this.eccentricity ** 2)) / (1.0 + this.eccentricity * Math.cos(this.trueAnomaly)));
        const vVec = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D(-Math.sin(this.trueAnomaly), this.eccentricity + Math.cos(this.trueAnomaly), 0.0);
        const vPQW = vVec.scale(Math.sqrt(this.mu / (this.semimajorAxis * (1 - this.eccentricity ** 2))));
        const position = rPQW
            .rotZ(-this.argPerigee)
            .rotX(-this.inclination)
            .rotZ(-this.rightAscension);
        const velocity = vPQW
            .rotZ(-this.argPerigee)
            .rotX(-this.inclination)
            .rotZ(-this.rightAscension);
        return { position, velocity };
    }
    /**
     * Converts the classical elements to equinoctial elements.
     * @returns The equinoctial elements.
     */
    toEquinoctialElements() {
        const I = this.inclination > Math.PI / 2 ? -1 : 1;
        const h = this.eccentricity * Math.sin(this.argPerigee + I * this.rightAscension);
        const k = this.eccentricity * Math.cos(this.argPerigee + I * this.rightAscension);
        const meanAnomaly = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.newtonNu)(this.eccentricity, this.trueAnomaly).m;
        const lambda = (meanAnomaly + this.argPerigee + I * this.rightAscension);
        const a = this.semimajorAxis;
        const p = Math.tan(0.5 * this.inclination) ** I * Math.sin(this.rightAscension);
        const q = Math.tan(0.5 * this.inclination) ** I * Math.cos(this.rightAscension);
        return new _EquinoctialElements_js__WEBPACK_IMPORTED_MODULE_5__.EquinoctialElements({ epoch: this.epoch, k, h, lambda, a, p, q, mu: this.mu, I });
    }
    /**
     * Propagates the classical elements to a given epoch.
     * @param propEpoch - The epoch to propagate the classical elements to.
     * @returns The classical elements at the propagated epoch.
     */
    propagate(propEpoch) {
        const t = this.epoch;
        const n = this.meanMotion;
        const delta = propEpoch.difference(t);
        const cosV = Math.cos(this.trueAnomaly);
        let eaInit = Math.acos((0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.clamp)((this.eccentricity + cosV) / (1 + this.eccentricity * cosV), -1, 1));
        eaInit = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.matchHalfPlane)(eaInit, this.trueAnomaly);
        let maInit = eaInit - this.eccentricity * Math.sin(eaInit);
        maInit = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.matchHalfPlane)(maInit, eaInit);
        const maFinal = (maInit + n * delta) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU;
        let eaFinal = maFinal;
        for (let iter = 0; iter < 32; iter++) {
            const eaTemp = maFinal + this.eccentricity * Math.sin(eaFinal);
            if (Math.abs(eaTemp - eaFinal) < 1e-12) {
                break;
            }
            eaFinal = eaTemp;
        }
        const cosEaFinal = Math.cos(eaFinal);
        let vFinal = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.clamp)(Math.acos((cosEaFinal - this.eccentricity) / (1 - this.eccentricity * cosEaFinal)), -1, 1);
        vFinal = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.matchHalfPlane)(vFinal, eaFinal);
        return new ClassicalElements({
            epoch: propEpoch,
            semimajorAxis: this.semimajorAxis,
            eccentricity: this.eccentricity,
            inclination: this.inclination,
            rightAscension: this.rightAscension,
            argPerigee: this.argPerigee,
            trueAnomaly: vFinal,
            mu: this.mu,
        });
    }
}


}),
"./src/engine/ootk/src/coordinate/EquinoctialElements.ts": 
/*!***************************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/EquinoctialElements.ts ***!
  \***************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EquinoctialElements: () => (EquinoctialElements)
});
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _ClassicalElements_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./ClassicalElements.js */ "./src/engine/ootk/src/coordinate/ClassicalElements.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



/**
 * Equinoctial elements are a set of orbital elements used to describe the
 * orbits of celestial bodies, such as satellites around a planet. They provide
 * an alternative to the traditional Keplerian elements and are especially
 * useful for avoiding singularities and numerical issues in certain types of
 * orbits.
 *
 * Unlike Keplerian elements, equinoctial elements don't suffer from
 * singularities at zero eccentricity (circular orbits) or zero inclination
 * (equatorial orbits). This makes them more reliable for numerical simulations
 * and analytical studies, especially in these edge cases.
 * @see https://faculty.nps.edu/dad/orbital/th0.pdf
 */
class EquinoctialElements {
    epoch;
    /** The semi-major axis of the orbit in kilometers. */
    a;
    /** The h component of the eccentricity vector. */
    h;
    /** The k component of the eccentricity vector. */
    k;
    /** The p component of the ascending node vector. */
    p;
    /** The q component of the ascending node vector. */
    q;
    /** The mean longitude of the orbit in radians. */
    lambda;
    /** The gravitational parameter of the central body in km³/s². */
    mu;
    /** The retrograde factor. 1 for prograde orbits, -1 for retrograde orbits. */
    I;
    constructor({ epoch, h, k, lambda, a, p, q, mu, I }) {
        this.epoch = epoch;
        this.h = h;
        this.k = k;
        this.lambda = lambda;
        this.a = a;
        this.p = p;
        this.q = q;
        this.mu = mu ?? _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.earthGravityParam;
        this.I = I ?? 1;
    }
    /**
     * Returns a string representation of the EquinoctialElements object.
     * @returns A string representation of the EquinoctialElements object.
     */
    toString() {
        return [
            '[EquinoctialElements]',
            `  Epoch: ${this.epoch}`,
            `  a: ${this.a} km`,
            `  h: ${this.h}`,
            `  k: ${this.k}`,
            `  p: ${this.p}`,
            `  q: ${this.q}`,
            `  lambda: ${this.lambda} rad`,
        ].join('\n');
    }
    /**
     * Gets the semimajor axis.
     * @returns The semimajor axis in kilometers.
     */
    get semimajorAxis() {
        return this.a;
    }
    /**
     * Gets the mean longitude.
     * @returns The mean longitude in radians.
     */
    get meanLongitude() {
        return this.lambda;
    }
    /**
     * Calculates the mean motion of the celestial object.
     * @returns The mean motion in units of radians per second.
     */
    get meanMotion() {
        return Math.sqrt(this.mu / this.a ** 3);
    }
    /**
     * Gets the retrograde factor.
     * @returns The retrograde factor.
     */
    get retrogradeFactor() {
        return this.I;
    }
    /**
     * Checks if the orbit is prograde.
     * @returns True if the orbit is prograde, false otherwise.
     */
    isPrograde() {
        return this.I === 1;
    }
    /**
     * Checks if the orbit is retrograde.
     * @returns True if the orbit is retrograde, false otherwise.
     */
    isRetrograde() {
        return this.I === -1;
    }
    /**
     * Gets the period of the orbit.
     * @returns The period in minutes.
     */
    get period() {
        const periodSec = (_utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.TAU * Math.sqrt(this.semimajorAxis ** 3 / this.mu));
        return (periodSec / 60);
    }
    /**
     * Gets the number of revolutions per day.
     * @returns The number of revolutions per day.
     */
    get revsPerDay() {
        return _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.MINUTES_PER_DAY / this.period;
    }
    /**
     * Converts the equinoctial elements to classical elements.
     * @returns The classical elements.
     */
    toClassicalElements() {
        const a = this.semimajorAxis;
        const e = Math.sqrt(this.k * this.k + this.h * this.h);
        const i = Math.PI * ((1.0 - this.I) * 0.5) + 2.0 * this.I * Math.atan(Math.sqrt(this.p * this.p + this.q * this.q));
        const o = Math.atan2(this.p, this.q);
        const w = Math.atan2(this.h, this.k) - this.I * Math.atan2(this.p, this.q);
        const m = this.lambda - this.I * o - w;
        const v = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_1__.newtonM)(e, m).nu;
        return new _ClassicalElements_js__WEBPACK_IMPORTED_MODULE_2__.ClassicalElements({
            epoch: this.epoch,
            semimajorAxis: a,
            eccentricity: e,
            inclination: i,
            rightAscension: o,
            argPerigee: w,
            trueAnomaly: v,
            mu: this.mu,
        });
    }
    /**
     * Converts the equinoctial elements to position and velocity.
     * @returns The position and velocity in classical elements.
     */
    toPositionVelocity() {
        return this.toClassicalElements().toPositionVelocity();
    }
}


}),
"./src/engine/ootk/src/coordinate/FormatTle.ts": 
/*!*****************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/FormatTle.ts ***!
  \*****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  FormatTle: () => (FormatTle)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * A class containing static methods for formatting TLEs (Two-Line Elements).
 */
class FormatTle {
    constructor() {
        // Static class
    }
    /**
     * Creates a TLE (Two-Line Element) string based on the provided TleParams.
     * @param tleParams - The parameters used to generate the TLE.
     * @returns An object containing the TLE strings tle1 and tle2.
     */
    static createTle(tleParams) {
        const { inc, meanmo, rasc, argPe, meana, ecen, epochyr, epochday, intl } = tleParams;
        const scc = _main_js__WEBPACK_IMPORTED_MODULE_0__.Tle.convert6DigitToA5(tleParams.scc);
        const epochYrStr = epochyr.padStart(2, '0');
        const epochdayStr = parseFloat(epochday).toFixed(8).padStart(12, '0');
        const incStr = FormatTle.inclination(inc);
        const meanmoStr = FormatTle.meanMotion(meanmo);
        const rascStr = FormatTle.rightAscension(rasc);
        const argPeStr = FormatTle.argumentOfPerigee(argPe);
        const meanaStr = FormatTle.meanAnomaly(meana);
        const ecenStr = FormatTle.eccentricity(ecen);
        const intlStr = intl.padEnd(8, ' ');
        // M' and M'' are both set to 0 to put the object in a perfect stable orbit
        let TLE1Ending = tleParams.sat ? tleParams.sat.tle1.substring(32, 71) : ' +.00000000 +00000+0 +00000-0 0  9990';
        // Add explicit positive/negative signs
        TLE1Ending = TLE1Ending[1] === ' ' ? FormatTle.setCharAt(TLE1Ending, 1, '+') : TLE1Ending;
        TLE1Ending = TLE1Ending[12] === ' ' ? FormatTle.setCharAt(TLE1Ending, 12, '+') : TLE1Ending;
        TLE1Ending = TLE1Ending[21] === ' ' ? FormatTle.setCharAt(TLE1Ending, 21, '+') : TLE1Ending;
        TLE1Ending = TLE1Ending[32] === ' ' ? FormatTle.setCharAt(TLE1Ending, 32, '0') : TLE1Ending;
        const tle1 = `1 ${scc}U ${intlStr} ${epochYrStr}${epochdayStr}${TLE1Ending}`;
        const tle2 = `2 ${scc} ${incStr} ${rascStr} ${ecenStr} ${argPeStr} ${meanaStr} ${meanmoStr} 00010`;
        return { tle1: tle1, tle2: tle2 };
    }
    /**
     * Converts the argument of perigee to a stringified number.
     * @param argPe - The argument of perigee to be converted. Can be either a number or a string.
     * @returns The argument of perigee as a stringified number.
     * @throws Error if the length of the argument of perigee is not 8.
     */
    static argumentOfPerigee(argPe) {
        if (typeof argPe === 'number') {
            argPe = argPe.toString();
        }
        const argPeNum = parseFloat(argPe).toFixed(4);
        const argPe0 = argPeNum.padStart(8, '0');
        if (argPe0.length !== 8) {
            throw new Error('argPe length is not 8');
        }
        return argPe0;
    }
    /**
     * Returns the eccentricity value of a given string.
     * @param ecen - The string representing the eccentricity.
     * @returns The eccentricity value.
     * @throws Error if the length of the eccentricity string is not 7.
     */
    static eccentricity(ecen) {
        let ecen0 = ecen.padEnd(9, '0');
        if (ecen0[1] === '.') {
            ecen0 = ecen0.substring(2);
        }
        else {
            ecen0 = ecen0.substring(0, 7);
        }
        if (ecen0.length !== 7) {
            throw new Error('ecen length is not 7');
        }
        return ecen0;
    }
    /**
     * Converts the inclination value to a string representation.
     * @param inc - The inclination value to be converted.
     * @returns The string representation of the inclination value.
     * @throws Error if the length of the converted value is not 8.
     */
    static inclination(inc) {
        if (typeof inc === 'number') {
            inc = inc.toString();
        }
        const incNum = parseFloat(inc).toFixed(4);
        const inc0 = incNum.padStart(8, '0');
        if (inc0.length !== 8) {
            throw new Error('inc length is not 8');
        }
        return inc0;
    }
    /**
     * Converts the mean anomaly to a string representation with 8 digits, padded with leading zeros.
     * @param meana - The mean anomaly to be converted. Can be either a number or a string.
     * @returns The mean anomaly as a string with 8 digits, padded with leading zeros.
     * @throws Error if the length of the mean anomaly is not 8.
     */
    static meanAnomaly(meana) {
        if (typeof meana === 'number') {
            meana = meana.toString();
        }
        const meanaNum = parseFloat(meana).toFixed(4);
        const meana0 = meanaNum.padStart(8, '0');
        if (meana0.length !== 8) {
            throw new Error('meana length is not 8');
        }
        return meana0;
    }
    /**
     * Converts the mean motion value to a string representation with 8 decimal
     * places. If the input is a number, it is converted to a string. If the input
     * is already a string, it is parsed as a float and then converted to a string
     * with 8 decimal places. The resulting string is padded with leading zeros to
     * ensure a length of 11 characters. Throws an error if the resulting string
     * does not have a length of 11 characters.
     * @param meanmo - The mean motion value to be converted.
     * @returns The string representation of the mean motion value with 8 decimal
     * places and padded with leading zeros.
     * @throws Error if the resulting string does not have a length of 11
     * characters.
     */
    static meanMotion(meanmo) {
        if (typeof meanmo === 'number') {
            meanmo = meanmo.toString();
        }
        const meanmoNum = parseFloat(meanmo).toFixed(8);
        const meanmo0 = meanmoNum.padStart(11, '0');
        if (meanmo0.length !== 11) {
            throw new Error('meanmo length is not 11');
        }
        return meanmo0;
    }
    /**
     * Converts the right ascension value to a stringified number.
     * @param rasc - The right ascension value to convert.
     * @returns The stringified number representation of the right ascension.
     * @throws Error if the length of the converted right ascension is not 8.
     */
    static rightAscension(rasc) {
        if (typeof rasc === 'number') {
            rasc = rasc.toString();
        }
        const rascNum = parseFloat(rasc).toFixed(4);
        const rasc0 = rascNum.padStart(8, '0');
        if (rasc0.length !== 8) {
            throw new Error('rasc length is not 8');
        }
        return rasc0;
    }
    /**
     * Sets a character at a specific index in a string. If the index is out of range, the original string is returned.
     * @param str - The input string.
     * @param index - The index at which to set the character.
     * @param chr - The character to set at the specified index.
     * @returns The modified string with the character set at the specified index.
     */
    static setCharAt(str, index, chr) {
        if (index > str.length - 1) {
            return str;
        }
        return `${str.substring(0, index)}${chr}${str.substring(index + 1)}`;
    }
}


}),
"./src/engine/ootk/src/coordinate/Geodetic.ts": 
/*!****************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/Geodetic.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Geodetic: () => (Geodetic)
});
/* ESM import */var _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../body/Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _ITRF_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./ITRF.js */ "./src/engine/ootk/src/coordinate/ITRF.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






/**
 * This Geodetic class represents a geodetic coordinate in three-dimensional
 * space, consisting of latitude, longitude, and altitude. It provides various
 * methods to perform calculations and operations related to geodetic
 * coordinates.
 *
 * This is a class for geodetic coordinates. This is related to the GroundObject
 * class, which is used to represent an object on the surface of the Earth.
 */
class Geodetic {
    lat;
    lon;
    alt;
    constructor(latitude, longitude, altitude) {
        if (Math.abs(latitude) > Math.PI / 2) {
            throw new RangeError('Latitude must be between -90° and 90° in Radians.');
        }
        if (Math.abs(longitude) > Math.PI) {
            throw new RangeError('Longitude must be between -180° and 180° in Radians.');
        }
        if (altitude < -_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean) {
            throw new RangeError(`Altitude must be greater than ${-_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean} km. Got ${altitude} km.`);
        }
        this.lat = latitude;
        this.lon = longitude;
        this.alt = altitude;
    }
    /**
     * Creates a Geodetic object from latitude, longitude, and altitude values in
     * degrees.
     * @param latitude The latitude value in degrees.
     * @param longitude The longitude value in degrees.
     * @param altitude The altitude value in kilometers.
     * @returns A Geodetic object representing the specified latitude, longitude,
     * and altitude.
     */
    static fromDegrees(latitude, longitude, altitude) {
        return new Geodetic((latitude * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD), (longitude * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD), altitude);
    }
    /**
     * Returns a string representation of the Geodetic object.
     * @returns A string containing the latitude, longitude, and altitude of the Geodetic object.
     */
    toString() {
        return [
            'Geodetic',
            `  Latitude:  ${this.latDeg.toFixed(4)}°`,
            `  Longitude: ${this.lonDeg.toFixed(4)}°`,
            `  Altitude:  ${this.alt.toFixed(3)} km`,
        ].join('\n');
    }
    /**
     * Gets the latitude in degrees.
     * @returns The latitude in degrees.
     */
    get latDeg() {
        return this.lat * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG;
    }
    /**
     * Gets the longitude in degrees.
     * @returns The longitude in degrees.
     */
    get lonDeg() {
        return this.lon * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG;
    }
    /**
     * Converts the geodetic coordinates to a ground position.
     * @returns The ground position object.
     */
    toGroundObject() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_1__.GroundObject({
            lat: this.latDeg,
            lon: this.lonDeg,
            alt: this.alt,
        });
    }
    /**
     * Converts the geodetic coordinates to the International Terrestrial
     * Reference Frame (ITRF) coordinates.
     * @param epoch The epoch in UTC.
     * @returns The ITRF coordinates.
     */
    toITRF(epoch) {
        const sLat = Math.sin(this.lat);
        const cLat = Math.cos(this.lat);
        const nVal = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator / Math.sqrt(1 - _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.eccentricitySquared * sLat * sLat);
        const r = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D(((nVal + this.alt) * cLat * Math.cos(this.lon)), ((nVal + this.alt) * cLat * Math.sin(this.lon)), ((nVal * (1 - _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.eccentricitySquared) + this.alt) * sLat));
        return new _ITRF_js__WEBPACK_IMPORTED_MODULE_5__.ITRF(epoch, r, _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D.origin);
    }
    /**
     * Calculates the angle between two geodetic coordinates.
     * @param g The geodetic coordinate to calculate the angle to.
     * @param method The method to use for calculating the angular distance (optional, default is Haversine).
     * @returns The angle between the two geodetic coordinates in radians.
     */
    angle(g, method = _main_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Haversine) {
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.angularDistance)(this.lon, this.lat, g.lon, g.lat, method);
    }
    /**
     * Calculates the angle in degrees between two Geodetic coordinates.
     * @param g The Geodetic coordinate to calculate the angle with.
     * @param method The method to use for calculating the angular distance (optional, default is Haversine).
     * @returns The angle in degrees.
     */
    angleDeg(g, method = _main_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Haversine) {
        return (this.angle(g, method) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG);
    }
    /**
     * Calculates the distance between two geodetic coordinates.
     * @param g The geodetic coordinates to calculate the distance to.
     * @param method The method to use for calculating the angular distance. Default is Haversine.
     * @returns The distance between the two geodetic coordinates in kilometers.
     */
    distance(g, method = _main_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Haversine) {
        return (this.angle(g, method) * _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean);
    }
    /**
     * Calculates the field of view based on the altitude of the Geodetic object.
     * @returns The field of view in radians.
     */
    fieldOfView() {
        return Math.acos(_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean / (_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean + this.alt));
    }
    /**
     * Determines if the current geodetic coordinate can see another geodetic coordinate.
     * @param g The geodetic coordinate to check for visibility.
     * @param method The method to use for calculating the angular distance (optional, default is Haversine).
     * @returns A boolean indicating if the current coordinate can see the other coordinate.
     */
    isInView(g, method = _main_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Haversine) {
        const fov = Math.max(this.fieldOfView(), g.fieldOfView());
        return this.angle(g, method) <= fov;
    }
}


}),
"./src/engine/ootk/src/coordinate/Hill.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/Hill.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Hill: () => (Hill)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _force_Thrust_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./../force/Thrust.js */ "./src/engine/ootk/src/force/Thrust.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


// / Hill Modified Equidistant Cyllindrical _(EQCM)_ coordinates.
class Hill {
    epoch;
    position;
    velocity;
    semimajorAxis_;
    meanMotion_;
    constructor(epoch, position, velocity, semimajorAxis) {
        this.epoch = epoch;
        this.position = position;
        this.velocity = velocity;
        this.semimajorAxis_ = semimajorAxis;
        this.meanMotion_ = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.smaToMeanMotion(this.semimajorAxis_);
    }
    static fromState(origin, radialPosition, intrackPosition, nodeVelocity, nodeOffsetTime) {
        const a = origin.semimajorAxis;
        const n = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.smaToMeanMotion(a);
        const yDot = -3.0 * radialPosition * n * 0.5;
        const z = (nodeVelocity / n) * Math.sin(n * -nodeOffsetTime);
        const zDot = nodeVelocity * Math.cos(n * -nodeOffsetTime);
        const r = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(radialPosition, intrackPosition, z);
        const v = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(0.0, yDot, zDot);
        return new Hill(origin.epoch, r, v, a);
    }
    static fromNmc(origin, majorAxisRange, nodeVelocity, nodeOffsetTime, translation = 0.0) {
        const a = origin.semimajorAxis;
        const n = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.smaToMeanMotion(a);
        const xDot = majorAxisRange * n * 0.5;
        const z = (nodeVelocity / n) * Math.sin(n * -nodeOffsetTime);
        const zDot = nodeVelocity * Math.cos(n * -nodeOffsetTime);
        const r = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(0.0, majorAxisRange + translation, z);
        const v = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(xDot, 0.0, zDot);
        return new Hill(origin.epoch, r, v, a);
    }
    static fromPerch(origin, perchRange, nodeVelocity, nodeOffsetTime) {
        const a = origin.semimajorAxis;
        const n = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.smaToMeanMotion(a);
        const z = (nodeVelocity / n) * Math.sin(n * -nodeOffsetTime);
        const zDot = nodeVelocity * Math.cos(n * -nodeOffsetTime);
        const r = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(0.0, perchRange, z);
        const v = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(0.0, 0.0, zDot);
        return new Hill(origin.epoch, r, v, a);
    }
    get semimajorAxis() {
        return this.semimajorAxis_;
    }
    set semimajorAxis(sma) {
        this.semimajorAxis_ = sma;
        this.meanMotion_ = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.smaToMeanMotion(this.semimajorAxis_);
    }
    get meanMotion() {
        return this.meanMotion_;
    }
    toJ2000Matrix(origin, transform) {
        const magrtgt = origin.position.magnitude();
        const magrint = magrtgt + this.position.x;
        const vtgtrsw = transform.multiplyVector3D(origin.velocity);
        const lambdadottgt = vtgtrsw.y / magrtgt;
        const lambdaint = this.position.y / magrtgt;
        const phiint = this.position.z / magrtgt;
        const sinphiint = Math.sin(phiint);
        const cosphiint = Math.cos(phiint);
        const sinlambdaint = Math.sin(lambdaint);
        const coslambdaint = Math.cos(lambdaint);
        const rotRswSez = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [sinphiint * coslambdaint, sinphiint * sinlambdaint, -cosphiint],
            [-sinlambdaint, coslambdaint, 0],
            [cosphiint * coslambdaint, cosphiint * sinlambdaint, sinphiint],
        ]);
        const rdotint = this.velocity.x + vtgtrsw.x;
        const lambdadotint = this.velocity.y / magrtgt + lambdadottgt;
        const phidotint = this.velocity.z / magrtgt;
        const vintsez = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(-magrint * phidotint, magrint * lambdadotint * cosphiint, rdotint);
        const vintrsw = rotRswSez.transpose().multiplyVector3D(vintsez);
        const vinteci = transform.transpose().multiplyVector3D(vintrsw);
        const rintrsw = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(cosphiint * magrint * coslambdaint, cosphiint * magrint * sinlambdaint, sinphiint * magrint);
        const rinteci = transform.transpose().multiplyVector3D(rintrsw);
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(origin.epoch, rinteci, vinteci);
    }
    toJ2000(origin) {
        return this.toJ2000Matrix(origin, _main_js__WEBPACK_IMPORTED_MODULE_0__.RelativeState.createMatrix(origin.position, origin.velocity));
    }
    static transitionMatrix(t, meanMotion) {
        const n = meanMotion;
        const sn = Math.sin(n * t);
        const cs = Math.cos(n * t);
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [4.0 - 3.0 * cs, 0.0, 0.0, sn / n, (2.0 * (1.0 - cs)) / n, 0.0],
            [6.0 * (sn - n * t), 1.0, 0.0, (-2.0 * (1.0 - cs)) / n, (4.0 * sn - 3.0 * n * t) / n, 0.0],
            [0.0, 0.0, cs, 0.0, 0.0, sn / n],
            [3.0 * n * sn, 0.0, 0.0, cs, 2.0 * sn, 0.0],
            [-6.0 * n * (1.0 - cs), 0.0, 0.0, -2.0 * sn, 4.0 * cs - 3.0, 0.0],
            [0.0, 0.0, -n * sn, 0.0, 0.0, cs],
        ]);
    }
    transition(t) {
        const sysMat = Hill.transitionMatrix(t, this.meanMotion_);
        const res = sysMat.multiplyVector(this.position.join(this.velocity)).elements;
        return new Hill(this.epoch.roll(t), new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(res[0], res[1], res[2]), new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(res[3], res[4], res[5]), this.semimajorAxis_);
    }
    transitionWithMatrix(stm, t) {
        const res = stm.multiplyVector(this.position.join(this.velocity)).elements;
        return new Hill(this.epoch.roll(t), new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(res[0], res[1], res[2]), new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(res[3], res[4], res[5]), this.semimajorAxis_);
    }
    propagate(newEpoch) {
        return this.transition(newEpoch.difference(this.epoch));
    }
    propagateWithMatrix(stm, newEpoch) {
        return this.transitionWithMatrix(stm, newEpoch.difference(this.epoch));
    }
    maneuver(maneuver) {
        const state = this.propagate(maneuver.center);
        return new Hill(state.epoch, state.position, state.velocity.add(maneuver.deltaV), state.semimajorAxis_);
    }
    ephemeris(start, stop, step = 60.0) {
        const output = [];
        let current = start;
        while (stop >= current) {
            output.push(this.propagate(current));
            current = current.roll(step);
        }
        return output;
    }
    get period() {
        return (2 * Math.PI) / this.meanMotion_;
    }
    nextRadialTangent() {
        const x = this.position.x;
        const xDot = this.velocity.x;
        const yDot = this.velocity.y;
        let t = Math.atan(-xDot / (3.0 * this.meanMotion_ * x + 2.0 * yDot)) / this.meanMotion_;
        if (t <= 0) {
            t = t + 0.5 * this.period;
        }
        else if (isNaN(t)) {
            t = 0.5 * this.period;
        }
        return this.propagate(this.epoch.roll(t));
    }
    solveManeuver(waypoint, ignoreCrosstrack = false) {
        const t = waypoint.epoch.difference(this.epoch);
        const w = waypoint.relativePosition;
        const sysMat = Hill.transitionMatrix(t, this.meanMotion_);
        const posEquationMat = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [sysMat.elements[0][0], sysMat.elements[0][1], sysMat.elements[0][2]],
            [sysMat.elements[1][0], sysMat.elements[1][1], sysMat.elements[1][2]],
            [sysMat.elements[2][0], sysMat.elements[2][1], sysMat.elements[2][2]],
        ]);
        const solnVector = w
            .subtract(posEquationMat.multiplyVector3D(this.position));
        const velEquationMat = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [sysMat.elements[0][3], sysMat.elements[0][4], sysMat.elements[0][5]],
            [sysMat.elements[1][3], sysMat.elements[1][4], sysMat.elements[1][5]],
            [sysMat.elements[2][3], sysMat.elements[2][4], sysMat.elements[2][5]],
        ]);
        let result = velEquationMat
            .inverse()
            .multiplyVector3D(solnVector)
            .subtract(this.velocity);
        if (ignoreCrosstrack) {
            result = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(result.x, result.y, 0);
        }
        return new _force_Thrust_js__WEBPACK_IMPORTED_MODULE_1__.Thrust(this.epoch, result.x * 1000, result.y * 1000, result.z * 1000);
    }
    maneuverSequence(pivot, waypoints, preManeuvers = [], postManeuvers = []) {
        let state = new Hill(this.epoch, this.position, this.velocity, this.semimajorAxis_);
        preManeuvers = preManeuvers.slice();
        postManeuvers = postManeuvers.slice();
        let output = preManeuvers;
        // Note difference was once compareTo
        output.sort((a, b) => a.center.difference(b.center));
        output = output.filter((mvr) => mvr.center >= this.epoch && mvr.center >= pivot);
        for (const mvr of output) {
            state = state.maneuver(mvr);
        }
        state = state.propagate(pivot);
        for (const wpt of waypoints) {
            const mvr = state.solveManeuver(wpt);
            state = state.maneuver(mvr);
            output.push(mvr);
        }
        output.push(...postManeuvers);
        return output;
    }
    maneuverOrigin(maneuver) {
        const state = this.propagate(maneuver.center);
        const vInit = Math.sqrt(_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / this.semimajorAxis_);
        const vFinal = vInit - maneuver.intrack * 1e-3;
        const aFinal = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / (vFinal * vFinal);
        return new Hill(state.epoch, state.position, state.velocity.subtract(maneuver.deltaV), aFinal);
    }
    get name() {
        return 'Hill';
    }
    toString() {
        return [
            `[${this.name}]`,
            `  Epoch: ${this.epoch}`,
            `  Position: ${this.position.toString(6)} km`,
            `  Velocity: ${this.velocity.toString(9)} km/s`,
        ].join('\n');
    }
}


}),
"./src/engine/ootk/src/coordinate/ITRF.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/ITRF.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ITRF: () => (ITRF)
});
/* ESM import */var _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../body/Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _Geodetic_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Geodetic.js */ "./src/engine/ootk/src/coordinate/Geodetic.ts");
/* ESM import */var _J2000_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./J2000.js */ "./src/engine/ootk/src/coordinate/J2000.ts");
/* ESM import */var _StateVector_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./StateVector.js */ "./src/engine/ootk/src/coordinate/StateVector.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable class-methods-use-this */




/**
 * The International Terrestrial Reference Frame (ITRF) is a geocentric reference frame for the Earth. It is the
 * successor to the International Terrestrial Reference System (ITRS). The ITRF definition is maintained by the
 * International Earth Rotation and Reference Systems Service (IERS). Several versions of ITRF exist, each with a
 * different epoch, to address the issue of crustal motion. The latest version is ITRF2014, based on data collected from
 * 1980 to 2014.
 * @see https://en.wikipedia.org/wiki/International_Terrestrial_Reference_Frame
 *
 * This is a geocentric coordinate system, also referenced as ECF/ECEF (Earth Centered Earth Fixed). It is a Cartesian
 * coordinate system with the origin at the center of the Earth. The x-axis intersects the sphere of the Earth at 0°
 * latitude (the equator) and 0° longitude (the Prime Meridian). The z-axis goes through the North Pole. The y-axis goes
 * through 90° East longitude.
 * @see https://en.wikipedia.org/wiki/Earth-centered,_Earth-fixed_coordinate_system
 */
class ITRF extends _StateVector_js__WEBPACK_IMPORTED_MODULE_3__.StateVector {
    /**
     * Gets the name of the ITRF coordinate system.
     * @returns The name of the coordinate system.
     */
    get name() {
        return 'ITRF';
    }
    /**
     * Gets a value indicating whether the coordinate system is inertial.
     * @returns A boolean value indicating whether the coordinate system is inertial.
     */
    get inertial() {
        return false;
    }
    /**
     * Gets the height of the ITRF coordinate above the surface of the Earth in kilometers.
     * @returns The height in kilometers.
     */
    get height() {
        const a = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator;
        const e2 = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.eccentricitySquared;
        const r = this.position.magnitude();
        const sl = this.position.z / r;
        const cl2 = 1 - sl * sl;
        const coeff = Math.sqrt((1 - e2) / (1 - e2 * cl2));
        return (r - a * coeff);
    }
    /**
     * Gets the altitude in kilometers.
     * @returns The altitude in kilometers.
     */
    get alt() {
        return this.height;
    }
    /**
     * Converts the current coordinate to the J2000 coordinate system. This is an Earth-Centered Inertial (ECI) coordinate
     * system with the origin at the center of the Earth.
     * @see https://en.wikipedia.org/wiki/Epoch_(astronomy)#Julian_years_and_J2000
     * @returns The coordinate in the J2000 coordinate system.
     */
    toJ2000() {
        const p = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.precession(this.epoch);
        const n = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.nutation(this.epoch);
        const ast = this.epoch.gmstAngle() + n.eqEq;
        const rTOD = this.position.rotZ(-ast);
        const vTOD = this.velocity
            .add(_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.rotation.cross(this.position))
            .rotZ(-ast);
        const rMOD = rTOD.rotX(n.eps).rotZ(n.dPsi).rotX(-n.mEps);
        const vMOD = vTOD.rotX(n.eps).rotZ(n.dPsi).rotX(-n.mEps);
        const rJ2000 = rMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
        const vJ2000 = vMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
        return new _J2000_js__WEBPACK_IMPORTED_MODULE_2__.J2000(this.epoch, rJ2000, vJ2000);
    }
    /**
     * Converts the current ITRF coordinate to Geodetic coordinate. This is a coordinate system for latitude, longitude,
     * and altitude.
     * @returns The converted Geodetic coordinate.
     */
    toGeodetic() {
        const sma = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator;
        const esq = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.eccentricitySquared;
        const x = this.position.x;
        const y = this.position.y;
        const z = this.position.z;
        const lon = Math.atan2(y, x);
        const r = Math.sqrt(x * x + y * y);
        const phi = Math.atan(z / r);
        let lat = phi;
        let alt;
        let c = 0.0;
        if (x === 0 && y === 0) {
            lat = phi;
            alt = z > 0 ? (z - _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusPolar) : (z + _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusPolar);
        }
        else {
            for (let i = 0; i < 20; i++) {
                const slat = Math.sin(lat);
                c = 1 / Math.sqrt(1 - esq * slat * slat);
                lat = Math.atan((z + sma * c * esq * slat) / r);
            }
            alt = (r / Math.cos(lat) - sma * c);
        }
        return new _Geodetic_js__WEBPACK_IMPORTED_MODULE_1__.Geodetic(lat, lon, alt);
    }
}


}),
"./src/engine/ootk/src/coordinate/J2000.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/J2000.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  J2000: () => (J2000)
});
/* ESM import */var _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../body/Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _ITRF_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./ITRF.js */ "./src/engine/ootk/src/coordinate/ITRF.ts");
/* ESM import */var _StateVector_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./StateVector.js */ "./src/engine/ootk/src/coordinate/StateVector.ts");
/* ESM import */var _TEME_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./TEME.js */ "./src/engine/ootk/src/coordinate/TEME.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




/**
 * Represents a position and velocity in the J2000 coordinate system. This is an Earth-centered inertial (ECI)
 * coordinate system.
 *
 * Commonly used ECI frame is defined with the Earth's Mean Equator and Mean Equinox (MEME) at 12:00 Terrestrial Time on
 * 1 January 2000. It can be referred to as J2K, J2000 or EME2000. The x-axis is aligned with the mean vernal equinox.
 * The z-axis is aligned with the Earth's rotation axis (or equivalently, the celestial North Pole) as it was at that
 * time. The y-axis is rotated by 90° East about the celestial equator.
 * @see https://en.wikipedia.org/wiki/Earth-centered_inertial
 */
class J2000 extends _StateVector_js__WEBPACK_IMPORTED_MODULE_2__.StateVector {
    /**
     * Creates a J2000 coordinate from classical elements.
     * @param elements The classical elements.
     * @returns The J2000 coordinate.
     */
    static fromClassicalElements(elements) {
        const rv = elements.toPositionVelocity();
        return new J2000(elements.epoch, rv.position, rv.velocity);
    }
    /**
     * Gets the name of the coordinate system.
     * @returns The name of the coordinate system.
     */
    get name() {
        return 'J2000';
    }
    /**
     * Gets a value indicating whether the coordinate system is inertial.
     * @returns A boolean value indicating whether the coordinate system is inertial.
     */
    get inertial() {
        return true;
    }
    /**
     * Converts the coordinates from J2000 to the International Terrestrial Reference Frame (ITRF).
     * This is an ECI to ECF transformation.
     * @returns The ITRF coordinates.
     */
    toITRF() {
        const p = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.precession(this.epoch);
        const n = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.nutation(this.epoch);
        const ast = (this.epoch.gmstAngle() + n.eqEq);
        const rMOD = this.position
            .rotZ(-p.zeta)
            .rotY(p.theta)
            .rotZ(-p.zed);
        const vMOD = this.velocity
            .rotZ(-p.zeta)
            .rotY(p.theta)
            .rotZ(-p.zed);
        const rTOD = rMOD
            .rotX(n.mEps)
            .rotZ(-n.dPsi)
            .rotX(-n.eps);
        const vTOD = vMOD
            .rotX(n.mEps)
            .rotZ(-n.dPsi)
            .rotX(-n.eps);
        const rPEF = rTOD.rotZ(ast);
        const vPEF = vTOD.rotZ(ast).add(_body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.rotation.negate().cross(rPEF));
        return new _ITRF_js__WEBPACK_IMPORTED_MODULE_1__.ITRF(this.epoch, rPEF, vPEF);
    }
    /**
     * Converts the J2000 coordinate to the TEME coordinate.
     * @returns The TEME coordinate.
     */
    toTEME() {
        const p = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.precession(this.epoch);
        const n = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.nutation(this.epoch);
        const eps = n.mEps + n.dEps;
        const dPsiCosEps = (n.dPsi * Math.cos(eps));
        const rMOD = this.position
            .rotZ(-p.zeta)
            .rotY(p.theta)
            .rotZ(-p.zed);
        const vMOD = this.velocity
            .rotZ(-p.zeta)
            .rotY(p.theta)
            .rotZ(-p.zed);
        const rTEME = rMOD
            .rotX(n.mEps)
            .rotZ(-n.dPsi)
            .rotX(-eps)
            .rotZ(dPsiCosEps);
        const vTEME = vMOD
            .rotX(n.mEps)
            .rotZ(-n.dPsi)
            .rotX(-eps)
            .rotZ(dPsiCosEps);
        return new _TEME_js__WEBPACK_IMPORTED_MODULE_3__.TEME(this.epoch, rTEME, vTEME);
    }
}


}),
"./src/engine/ootk/src/coordinate/RIC.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/coordinate/RIC.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RIC: () => (RIC)
});
/* ESM import */var _J2000_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./J2000.js */ "./src/engine/ootk/src/coordinate/J2000.ts");
/* ESM import */var _RelativeState_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./RelativeState.js */ "./src/engine/ootk/src/coordinate/RelativeState.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


/**
 * Represents a Radial-Intrack-Crosstrack (RIC) coordinates.
 */
class RIC extends _RelativeState_js__WEBPACK_IMPORTED_MODULE_1__.RelativeState {
    /**
     * Gets the name of the RIC coordinate system.
     * @returns The name of the RIC coordinate system.
     */
    get name() {
        return 'RIC';
    }
    /**
     * Creates a new RIC (Radial-Intrack-Crosstrack) coordinate from the J2000 state vectors.
     * @param state - The J2000 state vector.
     * @param origin - The J2000 state vector of the origin.
     * @param transform - The transformation matrix.
     * @returns The RIC coordinate.
     */
    static fromJ2000Matrix(state, origin, transform) {
        const dr = state.position.subtract(origin.position);
        const dv = state.velocity.subtract(origin.velocity);
        return new RIC(transform.multiplyVector3D(dr), transform.multiplyVector3D(dv));
    }
    /**
     * Creates a RIC (Radial-Intrack-Crosstrack) coordinate system from a J2000 state and origin.
     * @param state The J2000 state.
     * @param origin The J2000 origin.
     * @returns The RIC coordinate system.
     */
    static fromJ2000(state, origin) {
        return RIC.fromJ2000Matrix(state, origin, _RelativeState_js__WEBPACK_IMPORTED_MODULE_1__.RelativeState.createMatrix(origin.position, origin.velocity));
    }
    /**
     * Transforms the current RIC coordinate to the J2000 coordinate system using the provided origin and transform
     * matrix.
     * @param origin The origin J2000 coordinate.
     * @param transform The transformation matrix.
     * @returns The transformed J2000 coordinate.
     */
    toJ2000Matrix(origin, transform) {
        const tt = transform.transpose();
        const tr = tt.multiplyVector3D(this.position);
        const tv = tt.multiplyVector3D(this.velocity);
        return new _J2000_js__WEBPACK_IMPORTED_MODULE_0__.J2000(origin.epoch, origin.position.add(tr), origin.velocity.add(tv));
    }
    /**
     * Transforms the current RIC coordinate to the J2000 coordinate system using the provided origin.
     * @param origin The origin J2000 coordinate.
     * @returns The transformed J2000 coordinate.
     */
    toJ2000(origin) {
        return this.toJ2000Matrix(origin, _RelativeState_js__WEBPACK_IMPORTED_MODULE_1__.RelativeState.createMatrix(origin.position, origin.velocity));
    }
}


}),
"./src/engine/ootk/src/coordinate/RelativeState.ts": 
/*!*********************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/RelativeState.ts ***!
  \*********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RelativeState: () => (RelativeState)
});
/* ESM import */var _operations_Matrix_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../operations/Matrix.js */ "./src/engine/ootk/src/operations/Matrix.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Represents the relative state of an object in 3D space.
 */
class RelativeState {
    position;
    velocity;
    constructor(position, velocity) {
        this.position = position;
        this.velocity = velocity;
    }
    /**
     * Returns a string representation of the RelativeState object. The string includes the name, position, and velocity
     * of the object.
     * @returns A string representation of the RelativeState object.
     */
    toString() {
        return [
            `[${this.name}]`,
            `  Position: ${this.position.toString(6)} km`,
            `  Velocity: ${this.velocity.toString(9)} km/s`,
        ].join('\n');
    }
    /**
     * Creates a matrix based on the given position and velocity vectors. The matrix represents the relative state of an
     * object in 3D space.
     * @param position - The position vector.
     * @param velocity - The velocity vector.
     * @returns The matrix representing the relative state.
     */
    static createMatrix(position, velocity) {
        const ru = position.normalize();
        const cu = position.cross(velocity).normalize();
        const iu = cu.cross(ru).normalize();
        return new _operations_Matrix_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [ru.x, ru.y, ru.z],
            [iu.x, iu.y, iu.z],
            [cu.x, cu.y, cu.z],
        ]);
    }
    /**
     * Calculates the range of the relative state.
     * @returns The range in kilometers.
     */
    get range() {
        return this.position.magnitude();
    }
    /**
     * Calculates the range rate of the relative state. Range rate is the dot product of the position and velocity divided
     * by the range.
     * @returns The range rate in Kilometers per second.
     */
    get rangeRate() {
        return this.position.dot(this.velocity) / this.range;
    }
}


}),
"./src/engine/ootk/src/coordinate/StateVector.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/StateVector.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  StateVector: () => (StateVector)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * A state vector is a set of coordinates used to specify the position and
 * velocity of an object in a particular reference frame.
 */
class StateVector {
    epoch;
    position;
    velocity;
    constructor(epoch, position, velocity) {
        this.epoch = epoch;
        this.position = position;
        this.velocity = velocity;
    }
    /**
     * Returns a string representation of the StateVector object. The string includes the name, epoch, position, and
     * velocity.
     * @returns A string representation of the StateVector object.
     */
    toString() {
        return [
            `[${this.name}]`,
            `  Epoch: ${this.epoch}`,
            `  Position: ${this.position.toString(6)} km`,
            `  Velocity: ${this.velocity.toString(9)} km/s`,
        ].join('\n');
    }
    /**
     * Calculates the mechanical energy of the state vector.
     * @returns The mechanical energy value.
     */
    get mechanicalEnergy() {
        const r = this.position.magnitude();
        const v = this.velocity.magnitude();
        return v * v * 0.5 - _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / r;
    }
    /**
     * Calculates the semimajor axis of the state vector.
     * @returns The semimajor axis in kilometers.
     */
    get semimajorAxis() {
        const energy = this.mechanicalEnergy;
        return (-_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / (2.0 * energy));
    }
    /**
     * Gets the period of the state vector in minutes.
     * @returns The period in minutes.
     */
    get period() {
        const a = this.semimajorAxis;
        const periodSeconds = _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * Math.sqrt((a * a * a) / _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu);
        return (periodSeconds / 60.0);
    }
    /**
     * Gets the angular rate of the state vector.
     * @returns The angular rate.
     */
    get angularRate() {
        const a = this.semimajorAxis;
        return Math.sqrt(_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / (a * a * a));
    }
    /**
     * Converts the state vector to classical elements.
     * @param mu The gravitational parameter of the celestial body. Defaults to Earth's gravitational parameter.
     * @returns The classical elements corresponding to the state vector.
     * @throws Error if classical elements are undefined for fixed frames.
     */
    toClassicalElements(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        if (!this.inertial) {
            throw new Error('Classical elements are undefined for fixed frames.');
        }
        return _main_js__WEBPACK_IMPORTED_MODULE_0__.ClassicalElements.fromStateVector(this, mu);
    }
}


}),
"./src/engine/ootk/src/coordinate/TEME.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/TEME.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  TEME: () => (TEME)
});
/* ESM import */var _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../body/Earth.js */ "./src/engine/ootk/src/body/Earth.ts");
/* ESM import */var _J2000_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./J2000.js */ "./src/engine/ootk/src/coordinate/J2000.ts");
/* ESM import */var _StateVector_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./StateVector.js */ "./src/engine/ootk/src/coordinate/StateVector.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



/**
 * True Equator Mean Equinox (TEME) is a coordinate system commonly used in satellite tracking and orbit prediction. It
 * is a reference frame that defines the position and orientation of an object relative to the Earth's equator and
 * equinox.
 *
 * By using the True Equator Mean Equinox (TEME) coordinate system, we can accurately describe the position and motion
 * of satellites relative to the Earth's equator and equinox. This is particularly useful for tracking and predicting
 * satellite orbits in various applications, such as satellite communication, navigation, and remote sensing.
 */
class TEME extends _StateVector_js__WEBPACK_IMPORTED_MODULE_2__.StateVector {
    /**
     * Gets the name of the coordinate system.
     * @returns The name of the coordinate system.
     */
    get name() {
        return 'TEME';
    }
    /**
     * Gets a value indicating whether the coordinate is inertial.
     * @returns A boolean value indicating whether the coordinate is inertial.
     */
    get inertial() {
        return true;
    }
    /**
     * Creates a TEME (True Equator Mean Equinox) object from classical orbital elements.
     * @param elements - The classical orbital elements.
     * @returns A new TEME object.
     */
    static fromClassicalElements(elements) {
        const rv = elements.toPositionVelocity();
        return new TEME(elements.epoch, rv.position, rv.velocity);
    }
    /**
     * Converts the TEME (True Equator Mean Equinox) coordinates to J2000 coordinates.
     * @returns The J2000 coordinates.
     */
    toJ2000() {
        const p = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.precession(this.epoch);
        const n = _body_Earth_js__WEBPACK_IMPORTED_MODULE_0__.Earth.nutation(this.epoch);
        const eps = n.mEps + n.dEps;
        const dPsiCosEps = n.dPsi * Math.cos(eps);
        const rMOD = this.position
            .rotZ(-dPsiCosEps)
            .rotX(eps)
            .rotZ(n.dPsi)
            .rotX(-n.mEps);
        const vMOD = this.velocity
            .rotZ(-dPsiCosEps)
            .rotX(eps)
            .rotZ(n.dPsi)
            .rotX(-n.mEps);
        const rJ2K = rMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
        const vJ2K = vMOD
            .rotZ(p.zed)
            .rotY(-p.theta)
            .rotZ(p.zeta);
        return new _J2000_js__WEBPACK_IMPORTED_MODULE_1__.J2000(this.epoch, rJ2K, vJ2K);
    }
}


}),
"./src/engine/ootk/src/coordinate/Tle.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/coordinate/Tle.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Tle: () => (Tle)
});
/* ESM import */var _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/Sgp4OpsMode.js */ "./src/engine/ootk/src/enums/Sgp4OpsMode.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _sgp4_sgp4_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../sgp4/sgp4.js */ "./src/engine/ootk/src/sgp4/sgp4.ts");
/* ESM import */var _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../time/EpochUTC.js */ "./src/engine/ootk/src/time/EpochUTC.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _index_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./index.js */ "./src/engine/ootk/src/coordinate/index.ts");
/* ESM import */var _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./tle-format-data.js */ "./src/engine/ootk/src/coordinate/tle-format-data.ts");
/* eslint-disable max-lines */
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */








/**
 * Tle is a static class with a collection of methods for working with TLEs.
 */
class Tle {
    line1;
    line2;
    epoch;
    satnum;
    satrec_;
    /**
     * Mapping of alphabets to their corresponding numeric values.
     */
    static alpha5_ = {
        A: '10',
        B: '11',
        C: '12',
        D: '13',
        E: '14',
        F: '15',
        G: '16',
        H: '17',
        // I is skipped on purpose
        J: '18',
        K: '19',
        L: '20',
        M: '21',
        N: '22',
        // O is skipped on purpose
        P: '23',
        Q: '24',
        R: '25',
        S: '26',
        T: '27',
        U: '28',
        V: '29',
        W: '30',
        X: '31',
        Y: '32',
        Z: '33',
    };
    /** The argument of perigee field. */
    static argPerigee_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(35, 42);
    /** The BSTAR drag term field. */
    static bstar_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(54, 61);
    /** The checksum field. */
    static checksum_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(69, 69);
    /** The classification field. */
    static classification_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(8, 8);
    /** The eccentricity field. */
    static eccentricity_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(27, 33);
    /** The element set number field. */
    static elsetNum_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(65, 68);
    /** The ephemeris type field. */
    static ephemerisType_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(63, 63);
    /** The epoch day field. */
    static epochDay_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(21, 32);
    /** The epoch year field. */
    static epochYear_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(19, 20);
    /** The inclination field. */
    static inclination_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(9, 16);
    /** The international designator launch number field. */
    static intlDesLaunchNum_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(12, 14);
    /** The international designator launch piece field. */
    static intlDesLaunchPiece_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(15, 17);
    /** The international designator year field. */
    static intlDesYear_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(10, 11);
    /** The line number field. */
    static lineNumber_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(1, 1);
    /** The mean anomaly field. */
    static meanAnom_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(44, 51);
    /** The first derivative of the mean motion field. */
    static meanMoDev1_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(34, 43);
    /** The second derivative of the mean motion field. */
    static meanMoDev2_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(45, 52);
    /** The mean motion field. */
    static meanMo_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(53, 63);
    /** The right ascension of the ascending node field. */
    static rightAscension_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(18, 25);
    /** The revolution number field. */
    static revNum_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(64, 68);
    /** The satellite number field. */
    static satNum_ = new _tle_format_data_js__WEBPACK_IMPORTED_MODULE_7__.TleFormatData(3, 7);
    constructor(line1, line2, opsMode = _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.AFSPC, gravConst = _sgp4_sgp4_js__WEBPACK_IMPORTED_MODULE_2__.Sgp4GravConstants.wgs72) {
        this.line1 = line1;
        this.line2 = line2;
        this.epoch = Tle.parseEpoch_(line1.substring(18, 32));
        this.satnum = parseInt(Tle.convertA5to6Digit(line1.substring(2, 7)));
        this.satrec_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.Sgp4.createSatrec(line1, line2, gravConst, opsMode);
    }
    toString() {
        return `${this.line1}\n${this.line2}`;
    }
    /**
     * Gets the semimajor axis of the TLE.
     * @returns The semimajor axis value.
     */
    get semimajorAxis() {
        return Tle.tleSma_(this.line2);
    }
    /**
     * Gets the eccentricity of the TLE.
     * @returns The eccentricity value.
     */
    get eccentricity() {
        return Tle.tleEcc_(this.line2);
    }
    /**
     * Gets the inclination of the TLE.
     * @returns The inclination in degrees.
     */
    get inclination() {
        return Tle.tleInc_(this.line2);
    }
    /**
     * Gets the inclination in degrees.
     * @returns The inclination in degrees.
     */
    get inclinationDegrees() {
        return Tle.tleInc_(this.line2) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.RAD2DEG;
    }
    /**
     * Gets the apogee of the TLE (Two-Line Elements) object.
     * Apogee is the point in an orbit that is farthest from the Earth.
     * It is calculated as the product of the semimajor axis and (1 + eccentricity).
     * @returns The apogee value.
     */
    get apogee() {
        return this.semimajorAxis * (1 + this.eccentricity);
    }
    /**
     * Gets the perigee of the TLE (Two-Line Element Set).
     * The perigee is the point in the orbit of a satellite or other celestial body where it is closest to the Earth.
     * It is calculated as the product of the semimajor axis and the difference between 1 and the eccentricity.
     * @returns The perigee value.
     */
    get perigee() {
        return this.semimajorAxis * (1 - this.eccentricity);
    }
    /**
     * Gets the period of the TLE in minutes.
     * @returns The period of the TLE in minutes.
     */
    get period() {
        const periodSec = (_utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.TAU * Math.sqrt(this.semimajorAxis ** 3 / _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.earthGravityParam));
        return (periodSec / 60);
    }
    /**
     * Parses the epoch string and returns the corresponding EpochUTC object.
     * @param epochStr - The epoch string to parse.
     * @returns The parsed EpochUTC object.
     */
    static parseEpoch_(epochStr) {
        let year = parseInt(epochStr.substring(0, 2));
        if (year >= 57) {
            year += 1900;
        }
        else {
            year += 2000;
        }
        const days = parseFloat(epochStr.substring(2, 14)) - 1;
        return _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_3__.EpochUTC.fromDateTimeString(`${year}-01-01T00:00:00.000Z`).roll(days * _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.secondsPerDay);
    }
    static calcElsetAge(tle1, nowInput, outputUnits = 'days') {
        nowInput ??= new Date();
        const currentYearFull = nowInput.getUTCFullYear();
        const currentYearShort = currentYearFull % 100;
        const epochYearShort = parseInt(tle1.substring(18, 20), 10);
        const epochDayOfYear = parseFloat(tle1.substring(20, 32));
        let epochYearFull;
        if (epochYearShort <= currentYearShort) {
            epochYearFull = 2000 + epochYearShort;
        }
        else {
            epochYearFull = 1900 + epochYearShort;
        }
        const epochJday = epochDayOfYear + (epochYearFull * 365);
        const currentJday = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.getDayOfYear)() + (currentYearFull * 365);
        const currentTime = (nowInput.getUTCHours() * 3600 + nowInput.getUTCMinutes() * 60 +
            nowInput.getUTCSeconds()) / 86400;
        const daysOld = (currentJday + currentTime) - epochJday;
        switch (outputUnits) {
            case 'hours':
                return daysOld * 24;
            case 'minutes':
                return daysOld * 1440;
            case 'seconds':
                return daysOld * 86400;
            default:
                return daysOld;
        }
    }
    /**
     * Propagates the TLE (Two-Line Element Set) to a specific epoch and returns the TEME (True Equator Mean Equinox)
     * coordinates.
     * @param epoch The epoch to propagate the TLE to.
     * @returns The TEME coordinates at the specified epoch.
     * @throws Error if propagation fails.
     */
    propagate(epoch) {
        const r = new Float64Array(3);
        const v = new Float64Array(3);
        const stateVector = _main_js__WEBPACK_IMPORTED_MODULE_1__.Sgp4.propagate(this.satrec_, epoch.difference(this.epoch) / 60.0);
        if (!stateVector) {
            throw new Error('Propagation failed');
        }
        Tle.sv2rv_(stateVector, r, v);
        return new _index_js__WEBPACK_IMPORTED_MODULE_6__.TEME(epoch, new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(r[0], r[1], r[2]), new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(v[0], v[1], v[2]));
    }
    /**
     * Converts the state vector to position and velocity arrays.
     * @param stateVector - The state vector containing position and velocity information.
     * @param r - The array to store the position values.
     * @param v - The array to store the velocity values.
     */
    static sv2rv_(stateVector, r, v) {
        const pos = stateVector.position;
        const vel = stateVector.velocity;
        r[0] = pos.x;
        r[1] = pos.y;
        r[2] = pos.z;
        v[0] = vel.x;
        v[1] = vel.y;
        v[2] = vel.z;
    }
    /**
     * Returns the current state of the satellite in the TEME coordinate system.
     * @returns The current state of the satellite.
     */
    currentState_() {
        const r = new Float64Array(3);
        const v = new Float64Array(3);
        const stateVector = _main_js__WEBPACK_IMPORTED_MODULE_1__.Sgp4.propagate(this.satrec_, 0.0);
        Tle.sv2rv_(stateVector, r, v);
        return new _index_js__WEBPACK_IMPORTED_MODULE_6__.TEME(this.epoch, new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(r[0], r[1], r[2]), new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(v[0], v[1], v[2]));
    }
    /**
     * Gets the state of the TLE in the TEME coordinate system.
     * @returns The state of the TLE in the TEME coordinate system.
     */
    get state() {
        return this.currentState_();
    }
    /**
     * Calculates the Semi-Major Axis (SMA) from the second line of a TLE.
     * @param line2 The second line of the TLE.
     * @returns The Semi-Major Axis (SMA) in kilometers.
     */
    static tleSma_(line2) {
        const n = parseFloat(line2.substring(52, 63));
        return _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.earthGravityParam ** (1 / 3) / ((_utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.TAU * n) / _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.secondsPerDay) ** (2 / 3);
    }
    /**
     * Parses the eccentricity value from the second line of a TLE.
     * @param line2 The second line of the TLE.
     * @returns The eccentricity value.
     */
    static tleEcc_(line2) {
        return parseFloat(`0.${line2.substring(26, 33)}`);
    }
    /**
     * Calculates the inclination angle from the second line of a TLE.
     * @param line2 The second line of the TLE.
     * @returns The inclination angle in radians.
     */
    static tleInc_(line2) {
        return parseFloat(line2.substring(8, 16)) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.DEG2RAD;
    }
    /**
     * Creates a TLE (Two-Line Element) object from classical orbital elements.
     * @param elements - The classical orbital elements.
     * @returns A TLE object.
     */
    static fromClassicalElements(elements) {
        const { epochYr, epochDay } = elements.epoch.toEpochYearAndDay();
        const intl = '58001A  ';
        const scc = '00001';
        const tles = _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.createTle({
            inc: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.inclination(elements.inclinationDegrees),
            meanmo: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.meanMotion(elements.revsPerDay),
            ecen: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.eccentricity(elements.eccentricity.toFixed(7)),
            argPe: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.argumentOfPerigee(elements.argPerigeeDegrees),
            meana: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.meanAnomaly((0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.newtonNu)(elements.eccentricity, elements.trueAnomaly).m * _utils_constants_js__WEBPACK_IMPORTED_MODULE_4__.RAD2DEG),
            rasc: _index_js__WEBPACK_IMPORTED_MODULE_6__.FormatTle.rightAscension(elements.rightAscensionDegrees),
            epochday: epochDay,
            epochyr: epochYr,
            scc,
            intl,
        });
        return new Tle(tles.tle1, tles.tle2);
    }
    /**
     * Argument of perigee.
     * @see https://en.wikipedia.org/wiki/Argument_of_perigee
     * @example 69.9862
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The argument of perigee in degrees (0 to 360).
     */
    static argOfPerigee(tleLine2) {
        const argPe = parseFloat(tleLine2.substring(Tle.argPerigee_.start, Tle.argPerigee_.stop));
        if (!(argPe >= 0 && argPe <= 360)) {
            throw new Error(`Invalid argument of perigee: ${argPe}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(argPe, 4);
    }
    /**
     * BSTAR drag term (decimal point assumed).  Estimates the effects of atmospheric drag on the satellite's motion.
     * @see https://en.wikipedia.org/wiki/BSTAR
     * @example 0.000036771
     * @description ('36771-4' in the original Tle or 0.36771 * 10 ^ -4)
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The drag coefficient.
     */
    static bstar(tleLine1) {
        const BSTAR_PART_2 = Tle.bstar_.start + 1;
        const BSTAR_PART_3 = Tle.bstar_.start + 6;
        const BSTAR_PART_4 = Tle.bstar_.stop - 1;
        const bstarSymbol = tleLine1.substring(Tle.bstar_.start, BSTAR_PART_2);
        // Decimal place is assumed
        let bstar1 = parseFloat(`0.${tleLine1.substring(BSTAR_PART_2, BSTAR_PART_3)}`);
        const exponentSymbol = tleLine1.substring(BSTAR_PART_3, BSTAR_PART_4);
        let exponent = parseInt(tleLine1.substring(BSTAR_PART_4, Tle.bstar_.stop));
        if (exponentSymbol === '-') {
            exponent *= -1;
        }
        else if (exponentSymbol !== '+') {
            throw new Error(`Invalid BSTAR symbol: ${bstarSymbol}`);
        }
        bstar1 *= 10 ** exponent;
        if (bstarSymbol === '-') {
            bstar1 *= -1;
        }
        else if (bstarSymbol === '+' || bstarSymbol === ' ') {
            // Do nothing
        }
        else {
            throw new Error(`Invalid BSTAR symbol: ${bstarSymbol}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(bstar1, 14);
    }
    /**
     * Tle line 1 checksum (modulo 10), for verifying the integrity of this line of the Tle.
     * @example 3
     * @param tleLine The first line of the Tle to parse.
     * @returns The checksum value (0 to 9)
     */
    static checksum(tleLine) {
        return parseInt(tleLine.substring(Tle.checksum_.start, Tle.checksum_.stop));
    }
    /**
     * Returns the satellite classification.
     * Some websites like https://KeepTrack.space and Celestrak.org will embed
     * information in this field about the source of the Tle.
     * @example 'U'
     * unclassified
     * @example 'C'
     * confidential
     * @example 'S'
     * secret
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The satellite classification.
     */
    static classification(tleLine1) {
        return tleLine1.substring(Tle.classification_.start, Tle.classification_.stop);
    }
    /**
     * Orbital eccentricity, decimal point assumed. All artificial Earth satellites have an eccentricity between 0
     * (perfect circle) and 1 (parabolic orbit).
     * @example 0.0006317
     * (`0006317` in the original Tle)
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The eccentricity of the satellite (0 to 1)
     */
    static eccentricity(tleLine2) {
        const ecc = parseFloat(`0.${tleLine2.substring(Tle.eccentricity_.start, Tle.eccentricity_.stop)}`);
        if (!(ecc >= 0 && ecc <= 1)) {
            throw new Error(`Invalid eccentricity: ${ecc}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(ecc, 7);
    }
    /**
     * Tle element set number, incremented for each new Tle generated.
     * @see https://en.wikipedia.org/wiki/Two-line_element_set
     * @example 999
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The element number (1 to 999)
     */
    static elsetNum(tleLine1) {
        return parseInt(tleLine1.substring(Tle.elsetNum_.start, Tle.elsetNum_.stop));
    }
    /**
     * Private value - used by United States Space Force to reference the orbit model used to generate the Tle. Will
     * always be seen as zero externally (e.g. by "us", unless you are "them" - in which case, hello!).
     *
     * Starting in 2024, this may contain a 4 if the Tle was generated using the new SGP4-XP model. Until the source code
     * is released, there is no way to support that format in JavaScript or TypeScript.
     * @example 0
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The ephemeris type.
     */
    static ephemerisType(tleLine1) {
        const ephemerisType = parseInt(tleLine1.substring(Tle.ephemerisType_.start, Tle.ephemerisType_.stop));
        if (ephemerisType !== 0 && ephemerisType !== 4) {
            throw new Error('Invalid ephemeris type');
        }
        if (ephemerisType === 4) {
            throw new Error('SGP4-XP is not supported');
        }
        return ephemerisType;
    }
    /**
     * Fractional day of the year when the Tle was generated (Tle epoch).
     * @example 206.18396726
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The day of the year the Tle was generated. (1 to 365.99999999)
     */
    static epochDay(tleLine1) {
        const epochDay = parseFloat(tleLine1.substring(Tle.epochDay_.start, Tle.epochDay_.stop));
        if (epochDay < 1 || epochDay > 366.99999999) {
            throw new Error('Invalid epoch day');
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(epochDay, 8);
    }
    /**
     * Year when the Tle was generated (Tle epoch), last two digits.
     * @example 17
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The year the Tle was generated. (0 to 99)
     */
    static epochYear(tleLine1) {
        const epochYear = parseInt(tleLine1.substring(Tle.epochYear_.start, Tle.epochYear_.stop));
        if (epochYear < 0 || epochYear > 99) {
            throw new Error('Invalid epoch year');
        }
        return epochYear;
    }
    /**
     * Year when the Tle was generated (Tle epoch), four digits.
     * @example 2008
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The year the Tle was generated. (1957 to 2056)
     */
    static epochYearFull(tleLine1) {
        const epochYear = parseInt(tleLine1.substring(Tle.epochYear_.start, Tle.epochYear_.stop));
        if (epochYear < 0 || epochYear > 99) {
            throw new Error('Invalid epoch year');
        }
        if (epochYear < 57) {
            return epochYear + 2000;
        }
        return epochYear + 1900;
    }
    /**
     * Inclination relative to the Earth's equatorial plane in degrees. 0 to 90 degrees is a prograde orbit and 90 to 180
     * degrees is a retrograde orbit.
     * @example 51.6400
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The inclination of the satellite. (0 to 180)
     */
    static inclination(tleLine2) {
        const inc = parseFloat(tleLine2.substring(Tle.inclination_.start, Tle.inclination_.stop));
        if (inc < 0 || inc > 180) {
            throw new Error(`Invalid inclination: ${inc}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(inc, 4);
    }
    /**
     * International Designator (COSPAR ID)
     * @see https://en.wikipedia.org/wiki/International_Designator
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The International Designator.
     */
    static intlDes(tleLine1) {
        const year2 = this.intlDesYear(tleLine1);
        // Some TLEs don't have a year, so we can't generate an IntlDes
        if (isNaN(year2)) {
            return '';
        }
        const year4 = year2 < 57 ? year2 + 2000 : year2 + 1900;
        const launchNum = this.intlDesLaunchNum(tleLine1);
        const launchPiece = this.intlDesLaunchPiece(tleLine1);
        return `${year4}-${launchNum.toString().padStart(3, '0')}${launchPiece}`;
    }
    /**
     * International Designator (COSPAR ID): Launch number of the year.
     * @example 67
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The launch number of the International Designator. (1 to 999)
     */
    static intlDesLaunchNum(tleLine1) {
        return parseInt(tleLine1.substring(Tle.intlDesLaunchNum_.start, Tle.intlDesLaunchNum_.stop));
    }
    /**
     * International Designator  (COSPAR ID): Piece of the launch.
     * @example 'A'
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The launch piece of the International Designator. (A to ZZZ)
     */
    static intlDesLaunchPiece(tleLine1) {
        return tleLine1.substring(Tle.intlDesLaunchPiece_.start, Tle.intlDesLaunchPiece_.stop).trim();
    }
    /**
     * International Designator (COSPAR ID): Last 2 digits of launch year.
     * @example 98
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The year of the International Designator. (0 to 99)
     */
    static intlDesYear(tleLine1) {
        return parseInt(tleLine1.substring(Tle.intlDesYear_.start, Tle.intlDesYear_.stop));
    }
    /**
     * This should always return a 1 or a 2.
     * @example 1
     * @param tleLine The first line of the Tle to parse.
     * @returns The line number of the Tle. (1 or 2)
     */
    static lineNumber(tleLine) {
        const lineNum = parseInt(tleLine.substring(Tle.lineNumber_.start, Tle.lineNumber_.stop));
        if (lineNum !== 1 && lineNum !== 2) {
            throw new Error('Invalid line number');
        }
        return lineNum;
    }
    /**
     * Mean anomaly. Indicates where the satellite was located within its orbit at the time of the Tle epoch.
     * @see https://en.wikipedia.org/wiki/Mean_Anomaly
     * @example 25.2906
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The mean anomaly of the satellite. (0 to 360)
     */
    static meanAnomaly(tleLine2) {
        const meanA = parseFloat(tleLine2.substring(Tle.meanAnom_.start, Tle.meanAnom_.stop));
        if (!(meanA >= 0 && meanA <= 360)) {
            throw new Error(`Invalid mean anomaly: ${meanA}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(meanA, 4);
    }
    /**
     * First Time Derivative of the Mean Motion divided by two.  Defines how mean motion changes over time, so Tle
     * propagators can still be used to make reasonable guesses when times are distant from the original Tle epoch. This
     * is recorded in units of orbits per day per day.
     * @example 0.00001961
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The first derivative of the mean motion.
     */
    static meanMoDev1(tleLine1) {
        const meanMoDev1 = parseFloat(tleLine1.substring(Tle.meanMoDev1_.start, Tle.meanMoDev1_.stop));
        if (isNaN(meanMoDev1)) {
            throw new Error('Invalid first derivative of mean motion.');
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(meanMoDev1, 8);
    }
    /**
     * Second Time Derivative of Mean Motion divided by six (decimal point assumed). Measures rate of change in the Mean
     * Motion Dot so software can make reasonable guesses when times are distant from the original Tle epoch. Usually
     * zero, unless the satellite is manuevering or in a decaying orbit. This is recorded in units of orbits per day per
     * day per day.
     * @example 0
     * '00000-0' in the original Tle or 0.00000 * 10 ^ 0
     * @param tleLine1 The first line of the Tle to parse.
     * @returns The second derivative of the mean motion.
     */
    static meanMoDev2(tleLine1) {
        const meanMoDev2 = parseFloat(tleLine1.substring(Tle.meanMoDev2_.start, Tle.meanMoDev2_.stop));
        if (isNaN(meanMoDev2)) {
            throw new Error('Invalid second derivative of mean motion.');
        }
        // NOTE: Should this limit to a specific number of decimals?
        return meanMoDev2;
    }
    /**
     * Revolutions around the Earth per day (mean motion).
     * @see https://en.wikipedia.org/wiki/Mean_Motion
     * @example 15.54225995
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The mean motion of the satellite. (0 to 18)
     */
    static meanMotion(tleLine2) {
        const meanMo = parseFloat(tleLine2.substring(Tle.meanMo_.start, Tle.meanMo_.stop));
        if (!(meanMo > 0 && meanMo <= 18)) {
            throw new Error(`Invalid mean motion: ${meanMo}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(meanMo, 8);
    }
    /**
     * Calculates the period of a satellite orbit based on the given Tle line 2.
     * @example 92.53035747
     * @param tleLine2 The Tle line 2.
     * @returns The period of the satellite orbit in minutes.
     */
    static period(tleLine2) {
        const meanMo = Tle.meanMotion(tleLine2);
        return (1440 / meanMo);
    }
    /**
     * Right ascension of the ascending node in degrees. Essentially, this is the angle of the satellite as it crosses
     * northward (ascending) across the Earth's equator (equatorial plane).
     * @example 208.9163
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The right ascension of the satellite. (0 to 360)
     */
    static rightAscension(tleLine2) {
        const rightAscension = parseFloat(tleLine2.substring(Tle.rightAscension_.start, Tle.rightAscension_.stop));
        if (!(rightAscension >= 0 && rightAscension <= 360)) {
            throw new Error(`Invalid Right Ascension: ${rightAscension}`);
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision)(rightAscension, 4);
    }
    /**
     * NORAD catalog number. To support Alpha-5, the first digit can be a letter. This will NOT be converted to a number.
     * Use satNum() for that.
     * @see https://en.wikipedia.org/wiki/Satellite_Catalog_Number
     * @example 25544
     * @example B1234
     * @param tleLine The first line of the Tle to parse.
     * @returns NORAD catalog number.
     */
    static rawSatNum(tleLine) {
        return tleLine.substring(Tle.satNum_.start, Tle.satNum_.stop);
    }
    /**
     * Total satellite revolutions when this Tle was generated. This number rolls over (e.g. 99999 -> 0).
     * @example 6766
     * @param tleLine2 The second line of the Tle to parse.
     * @returns The revolutions around the Earth per day (mean motion). (0 to 99999)
     */
    static revNum(tleLine2) {
        return parseInt(tleLine2.substring(Tle.revNum_.start, Tle.revNum_.stop));
    }
    /**
     * NORAD catalog number converted to a number.
     * @see https://en.wikipedia.org/wiki/Satellite_Catalog_Number
     * @example 25544
     * @example 111234
     * @param tleLine The first line of the Tle to parse.
     * @returns NORAD catalog number. (0 to 339999)
     */
    static satNum(tleLine) {
        const satNumStr = tleLine.substring(Tle.satNum_.start, Tle.satNum_.stop);
        const sixDigitSatNum = Tle.convertA5to6Digit(satNumStr);
        return parseInt(sixDigitSatNum);
    }
    /**
     * Parse the first line of the Tle.
     * @param tleLine1 The first line of the Tle to parse.
     * @returns Returns the data from the first line of the Tle.
     */
    static parseLine1(tleLine1) {
        const lineNumber1 = Tle.lineNumber(tleLine1);
        const satNum = Tle.satNum(tleLine1);
        const satNumRaw = Tle.rawSatNum(tleLine1);
        const classification = Tle.classification(tleLine1);
        const intlDes = Tle.intlDes(tleLine1);
        const intlDesYear = Tle.intlDesYear(tleLine1);
        const intlDesLaunchNum = Tle.intlDesLaunchNum(tleLine1);
        const intlDesLaunchPiece = Tle.intlDesLaunchPiece(tleLine1);
        const epochYear = Tle.epochYear(tleLine1);
        const epochYearFull = Tle.epochYearFull(tleLine1);
        const epochDay = Tle.epochDay(tleLine1);
        const meanMoDev1 = Tle.meanMoDev1(tleLine1);
        const meanMoDev2 = Tle.meanMoDev2(tleLine1);
        const bstar = Tle.bstar(tleLine1);
        const ephemerisType = Tle.ephemerisType(tleLine1);
        const elsetNum = Tle.elsetNum(tleLine1);
        const checksum1 = Tle.checksum(tleLine1);
        return {
            lineNumber1,
            satNum,
            satNumRaw,
            classification,
            intlDes,
            intlDesYear,
            intlDesLaunchNum,
            intlDesLaunchPiece,
            epochYear,
            epochYearFull,
            epochDay,
            meanMoDev1,
            meanMoDev2,
            bstar,
            ephemerisType,
            elsetNum,
            checksum1,
        };
    }
    /**
     * Parse the second line of the Tle.
     * @param tleLine2 The second line of the Tle to parse.
     * @returns Returns the data from the second line of the Tle.
     */
    static parseLine2(tleLine2) {
        const lineNumber2 = Tle.lineNumber(tleLine2);
        const satNum = Tle.satNum(tleLine2);
        const satNumRaw = Tle.rawSatNum(tleLine2);
        const inclination = Tle.inclination(tleLine2);
        const rightAscension = Tle.rightAscension(tleLine2);
        const eccentricity = Tle.eccentricity(tleLine2);
        const argOfPerigee = Tle.argOfPerigee(tleLine2);
        const meanAnomaly = Tle.meanAnomaly(tleLine2);
        const meanMotion = Tle.meanMotion(tleLine2);
        const revNum = Tle.revNum(tleLine2);
        const checksum2 = Tle.checksum(tleLine2);
        const period = Tle.period(tleLine2);
        return {
            lineNumber2,
            satNum,
            satNumRaw,
            inclination,
            rightAscension,
            eccentricity,
            argOfPerigee,
            meanAnomaly,
            meanMotion,
            revNum,
            checksum2,
            period,
        };
    }
    /**
     * Parses the Tle into orbital data.
     *
     * If you want all of the data then use parseTleFull instead.
     * @param tleLine1 Tle line 1
     * @param tleLine2 Tle line 2
     * @returns Returns most commonly used orbital data from Tle
     */
    static parse(tleLine1, tleLine2) {
        const line1 = Tle.parseLine1(tleLine1);
        const line2 = Tle.parseLine2(tleLine2);
        if (line1.satNum !== line2.satNum) {
            throw new Error('Satellite numbers do not match');
        }
        if (line1.satNumRaw !== line2.satNumRaw) {
            throw new Error('Raw satellite numbers do not match');
        }
        if (line1.lineNumber1 !== 1) {
            throw new Error('First line number must be 1');
        }
        if (line2.lineNumber2 !== 2) {
            throw new Error('Second line number must be 2');
        }
        return {
            satNum: line1.satNum,
            intlDes: line1.intlDes,
            epochYear: line1.epochYear,
            epochDay: line1.epochDay,
            meanMoDev1: line1.meanMoDev1,
            meanMoDev2: line1.meanMoDev2,
            bstar: line1.bstar,
            inclination: line2.inclination,
            rightAscension: line2.rightAscension,
            eccentricity: line2.eccentricity,
            argOfPerigee: line2.argOfPerigee,
            meanAnomaly: line2.meanAnomaly,
            meanMotion: line2.meanMotion,
            period: line2.period,
        };
    }
    /**
     * Parses all of the data contained in the Tle.
     *
     * If you only want the most commonly used data then use parseTle instead.
     * @param tleLine1 The first line of the Tle to parse.
     * @param tleLine2 The second line of the Tle to parse.
     * @returns Returns all of the data from the Tle.
     */
    static parseAll(tleLine1, tleLine2) {
        const line1 = Tle.parseLine1(tleLine1);
        const line2 = Tle.parseLine2(tleLine2);
        if (line1.satNum !== line2.satNum) {
            throw new Error('Satellite numbers do not match');
        }
        if (line1.satNumRaw !== line2.satNumRaw) {
            throw new Error('Raw satellite numbers do not match');
        }
        if (line1.lineNumber1 !== 1) {
            throw new Error('First line number must be 1');
        }
        if (line2.lineNumber2 !== 2) {
            throw new Error('Second line number must be 2');
        }
        return { ...line1, ...line2 };
    }
    /**
     * Converts a 6 digit SCC number to a 5 digit SCC alpha 5 number
     * @param sccNum The 6 digit SCC number
     * @returns The 5 digit SCC alpha 5 number
     */
    static convert6DigitToA5(sccNum) {
        // Only applies to 6 digit numbers
        if (sccNum.length < 6) {
            return sccNum;
        }
        if (typeof sccNum[0] !== 'string') {
            throw new Error('Invalid SCC number');
        }
        // Already an alpha 5 number
        if (RegExp(/[A-Z]/iu, 'u').test(sccNum[0])) {
            return sccNum;
        }
        // Extract the trailing 4 digits
        const rest = sccNum.slice(2, 6);
        /*
         * Convert the first two digit numbers into a Letter. Skip I and O as they
         * look too similar to 1 and 0 A=10, B=11, C=12, D=13, E=14, F=15, G=16,
         * H=17, J=18, K=19, L=20, M=21, N=22, P=23, Q=24, R=25, S=26, T=27, U=28,
         * V=29, W=30, X=31, Y=32, Z=33
         */
        let first = parseInt(`${sccNum[0]}${sccNum[1]}`);
        const iPlus = first >= 18 ? 1 : 0;
        const tPlus = first >= 24 ? 1 : 0;
        first = first + iPlus + tPlus;
        return `${String.fromCharCode(first + 55)}${rest}`;
    }
    /**
     * Converts a 5-digit SCC number to a 6-digit SCC number.
     * @param sccNum - The 5-digit SCC number to convert.
     * @returns The converted 6-digit SCC number.
     */
    static convertA5to6Digit(sccNum) {
        if (sccNum.length < 5) {
            return sccNum;
        }
        const values = sccNum.toUpperCase().split('');
        if (!values[0]) {
            throw new Error('Invalid SCC number');
        }
        if (values[0] in Tle.alpha5_) {
            const firstLetter = values[0];
            values[0] = Tle.alpha5_[firstLetter];
        }
        return values.join('');
    }
}


}),
"./src/engine/ootk/src/coordinate/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ClassicalElements: () => (/* reexport safe */ _ClassicalElements_js__WEBPACK_IMPORTED_MODULE_0__.ClassicalElements),
  EquinoctialElements: () => (/* reexport safe */ _EquinoctialElements_js__WEBPACK_IMPORTED_MODULE_1__.EquinoctialElements),
  FormatTle: () => (/* reexport safe */ _FormatTle_js__WEBPACK_IMPORTED_MODULE_2__.FormatTle),
  Geodetic: () => (/* reexport safe */ _Geodetic_js__WEBPACK_IMPORTED_MODULE_3__.Geodetic),
  Hill: () => (/* reexport safe */ _Hill_js__WEBPACK_IMPORTED_MODULE_11__.Hill),
  ITRF: () => (/* reexport safe */ _ITRF_js__WEBPACK_IMPORTED_MODULE_4__.ITRF),
  J2000: () => (/* reexport safe */ _J2000_js__WEBPACK_IMPORTED_MODULE_5__.J2000),
  RIC: () => (/* reexport safe */ _RIC_js__WEBPACK_IMPORTED_MODULE_7__.RIC),
  RelativeState: () => (/* reexport safe */ _RelativeState_js__WEBPACK_IMPORTED_MODULE_6__.RelativeState),
  StateVector: () => (/* reexport safe */ _StateVector_js__WEBPACK_IMPORTED_MODULE_8__.StateVector),
  TEME: () => (/* reexport safe */ _TEME_js__WEBPACK_IMPORTED_MODULE_9__.TEME),
  Tle: () => (/* reexport safe */ _Tle_js__WEBPACK_IMPORTED_MODULE_10__.Tle)
});
/* ESM import */var _ClassicalElements_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./ClassicalElements.js */ "./src/engine/ootk/src/coordinate/ClassicalElements.ts");
/* ESM import */var _EquinoctialElements_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./EquinoctialElements.js */ "./src/engine/ootk/src/coordinate/EquinoctialElements.ts");
/* ESM import */var _FormatTle_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./FormatTle.js */ "./src/engine/ootk/src/coordinate/FormatTle.ts");
/* ESM import */var _Geodetic_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Geodetic.js */ "./src/engine/ootk/src/coordinate/Geodetic.ts");
/* ESM import */var _ITRF_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./ITRF.js */ "./src/engine/ootk/src/coordinate/ITRF.ts");
/* ESM import */var _J2000_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./J2000.js */ "./src/engine/ootk/src/coordinate/J2000.ts");
/* ESM import */var _RelativeState_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./RelativeState.js */ "./src/engine/ootk/src/coordinate/RelativeState.ts");
/* ESM import */var _RIC_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./RIC.js */ "./src/engine/ootk/src/coordinate/RIC.ts");
/* ESM import */var _StateVector_js__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./StateVector.js */ "./src/engine/ootk/src/coordinate/StateVector.ts");
/* ESM import */var _TEME_js__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./TEME.js */ "./src/engine/ootk/src/coordinate/TEME.ts");
/* ESM import */var _Tle_js__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./Tle.js */ "./src/engine/ootk/src/coordinate/Tle.ts");
/* ESM import */var _Hill_js__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(/*! ./Hill.js */ "./src/engine/ootk/src/coordinate/Hill.ts");














}),
"./src/engine/ootk/src/coordinate/tle-format-data.ts": 
/*!***********************************************************!*\
  !*** ./src/engine/ootk/src/coordinate/tle-format-data.ts ***!
  \***********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  TleFormatData: () => (TleFormatData)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/**
 * Represents the format data of a TLE (Two-Line Element) set. This is used to
 * make it easier to remember the starting and ending positions of the columns
 * containing the TLE data.
 */
class TleFormatData {
    /** The starting position of the TLE data in the source string. */
    start;
    /** The ending position of the TLE data in the source string. */
    stop;
    /** The length of the TLE data in the source string. */
    length;
    /**
     * Creates a new instance of TleFormatData.
     * @param start The starting position of the TLE data in the source string.
     * @param end The ending position of the TLE data in the source string.
     */
    constructor(start, end) {
        this.start = start - 1;
        this.stop = end;
        this.length = this.stop - this.start;
    }
}


}),
"./src/engine/ootk/src/covariance/CovarianceSample.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/covariance/CovarianceSample.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CovarianceSample: () => (CovarianceSample)
});
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../propagator/RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/* ESM import */var _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./StateCovariance.js */ "./src/engine/ootk/src/covariance/StateCovariance.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable class-methods-use-this */




// / Sigma point covariance sample.
class CovarianceSample {
    origin_;
    samples_ = [];
    matrix_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.Matrix.zero(6, 12);
    /**
     * Create a new [CovarianceSample] object from an inertial state, covariance
     * and optional force models for the origin state and samples.
     *
     * Two-body physics will be used if a force model is not provided.
     * @param state The origin state.
     * @param covariance The covariance.
     * @param tle The TLE object.
     * @param originForceModel The force model for the origin state.
     * @param sampleForceModel The force model for the samples.
     */
    constructor(state, covariance, tle, originForceModel, sampleForceModel) {
        originForceModel ??= new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity();
        sampleForceModel ??= new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity();
        // Scale covariance using TLE quality and regime aging factor if TLE info is provided
        let scale = [1, 1, 1];
        if (tle) {
            const tleAgeDays = _main_js__WEBPACK_IMPORTED_MODULE_1__.Tle.calcElsetAge(tle.line1, new Date(), 'days');
            const quality = this.evaluateTleQuality(tle);
            const aging = this.getRegimeAgingFactor(tle, tleAgeDays);
            scale = [
                quality[0] * aging[0],
                quality[1] * aging[1],
                quality[2] * aging[2],
            ];
        }
        this.origin_ = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(state, originForceModel);
        const s = covariance.matrix.cholesky().elements;
        const sqrt6 = Math.sqrt(6.0);
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 6; j++) {
                /*
                 * Apply scale[0] to R, scale[1] to I, scale[2] to C (x, y, z)
                 * Position: i = 0,1,2; Velocity: i = 3,4,5
                 */
                const scaleIdx = i % 3;
                s[i][j] *= sqrt6 * scale[scaleIdx];
            }
        }
        // 6 x 12 matrix
        const sigmapts = _main_js__WEBPACK_IMPORTED_MODULE_1__.Matrix.zero(6, 12).elements;
        for (let i = 0; i < 6; i++) {
            const jj = (i - 1) * 2 + 2;
            for (let j = 0; j < 3; j++) {
                sigmapts[j][jj] = s[j][i];
                sigmapts[j + 3][jj] = s[j + 3][i];
                sigmapts[j][jj + 1] = -s[j][i];
                sigmapts[j + 3][jj + 1] = -s[j + 3][i];
            }
        }
        for (let i = 0; i < 12; i++) {
            const sampleR = new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(sigmapts[0][i], sigmapts[1][i], sigmapts[2][i]);
            const sampleV = new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(sigmapts[3][i], sigmapts[4][i], sigmapts[5][i]);
            if (covariance.frame === _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.CovarianceFrame.ECI) {
                const sample = new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(state.epoch, state.position.add(sampleR), state.velocity.add(sampleV));
                this.samples_.push(new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(sample, sampleForceModel));
            }
            else if (covariance.frame === _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.CovarianceFrame.RIC) {
                const sample = new _main_js__WEBPACK_IMPORTED_MODULE_1__.RIC(sampleR, sampleV).toJ2000(state);
                this.samples_.push(new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(sample, sampleForceModel));
            }
        }
    }
    // / Current covariance sample epoch.
    get epoch() {
        return this.origin_.state.epoch;
    }
    // / Current covariance sample origin state.
    get state() {
        return this.origin_.state;
    }
    // / Rebuild covariance from sigma points.
    _rebuildCovariance(matrix) {
        const pts = matrix.elements;
        const c = 1.0 / 12.0;
        const yu = new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector([0, 0, 0, 0, 0, 0]).elements;
        const y = _main_js__WEBPACK_IMPORTED_MODULE_1__.Matrix.zero(6, 12);
        for (let i = 0; i < 12; i++) {
            for (let j = 0; j < 6; j++) {
                yu[j] += pts[j][i];
            }
        }
        for (let j = 0; j < 6; j++) {
            yu[j] *= c;
        }
        for (let i = 0; i < 12; i++) {
            for (let j = 0; j < 6; j++) {
                y.elements[j][i] = pts[j][i] - yu[j];
            }
        }
        const yt = y.transpose();
        const tmp = y.multiply(yt);
        return tmp.scale(c);
    }
    // / Propagate covariance to a new epoch.
    propagate(epoch) {
        this.origin_.propagate(epoch);
        for (const sample of this.samples_) {
            sample.propagate(epoch);
        }
    }
    // / Apply a maneuver to this covariance.
    maneuver(maneuver) {
        this.origin_.maneuver(maneuver);
        for (const sample of this.samples_) {
            sample.maneuver(maneuver);
        }
    }
    // / Desample covariance in J2000 frame.
    desampleJ2000() {
        for (let i = 0; i < 12; i++) {
            const state = this.samples_[i].state;
            this.matrix_.elements[0][i] = state.position.x;
            this.matrix_.elements[1][i] = state.position.y;
            this.matrix_.elements[2][i] = state.position.z;
            this.matrix_.elements[3][i] = state.velocity.x;
            this.matrix_.elements[4][i] = state.velocity.y;
            this.matrix_.elements[5][i] = state.velocity.z;
        }
        const matrix = this._rebuildCovariance(this.matrix_);
        return new _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.StateCovariance(matrix, _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.CovarianceFrame.ECI);
    }
    // / Desample covariance in RIC frame.
    desampleRIC() {
        const rot = _main_js__WEBPACK_IMPORTED_MODULE_1__.RelativeState.createMatrix(this.origin_.state.position, this.origin_.state.velocity);
        for (let i = 0; i < 12; i++) {
            const state = _main_js__WEBPACK_IMPORTED_MODULE_1__.RIC.fromJ2000Matrix(this.samples_[i].state, this.origin_.state, rot);
            this.matrix_.elements[0][i] = state.position.x;
            this.matrix_.elements[1][i] = state.position.y;
            this.matrix_.elements[2][i] = state.position.z;
            this.matrix_.elements[3][i] = state.velocity.x;
            this.matrix_.elements[4][i] = state.velocity.y;
            this.matrix_.elements[5][i] = state.velocity.z;
        }
        const matrix = this._rebuildCovariance(this.matrix_);
        return new _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.StateCovariance(matrix, _StateCovariance_js__WEBPACK_IMPORTED_MODULE_3__.CovarianceFrame.RIC);
    }
    evaluateTleQuality(tle) {
        let c = 1, i = 1, r = 1; // start with nominal 120 / 1000 / 100 m
        /* ---- Mean–motion first derivative (rev/day²) ---- */
        const mm1 = Math.abs(_main_js__WEBPACK_IMPORTED_MODULE_1__.Tle.meanMoDev1(tle.line1));
        if (mm1 > 1e-3) { // clear manoeuvre or very low-drag orbit
            i *= 3.0;
            r *= 2.0;
            c *= 1.2;
        }
        else if (mm1 > 1e-5) { // high drag but likely passive
            i *= 1.6;
            r *= 1.3;
        }
        /* ---- BSTAR drag term ---- */
        const bstar = Math.abs(_main_js__WEBPACK_IMPORTED_MODULE_1__.Tle.bstar(tle.line1));
        if (bstar > 5e-4) { // <≈400 km LEO
            i *= 1.5;
            r *= 1.5;
            c *= 1.1;
        }
        else if (bstar > 1e-4) { // 400–600 km
            i *= 1.2;
            r *= 1.2;
        }
        /* ---- Eccentricity ---- */
        const e = tle.eccentricity;
        if (e > 0.05) { // HEO, GTO, cubesat transfer, etc.
            r *= 1 + 6 * e;
            i *= 1 + 4 * e;
            c *= 1 + 1 * e;
        }
        else if (e > 0.02) { // mildly elliptical LEO/MEO
            r *= 1 + 3 * e;
            i *= 1 + 2 * e;
        }
        return [r, i, c];
    }
    getRegimeAgingFactor(tle, ageDays) {
        const regime = tle.state.toClassicalElements().getOrbitRegime();
        const t = Math.max(ageDays, 0) ** 1.5; // non-linear growth
        switch (regime) {
            case _main_js__WEBPACK_IMPORTED_MODULE_1__.OrbitRegime.LEO:
                // Target ~×2 (R), ×3 (I), ×2.2 (C) after 1 day
                return [1 + 1.0 * t, 1 + 2.0 * t, 1 + 1.2 * t];
            case _main_js__WEBPACK_IMPORTED_MODULE_1__.OrbitRegime.MEO:
                // GNSS shells – slower growth
                return [1 + 0.4 * t, 1 + 0.9 * t, 1 + 0.6 * t];
            case _main_js__WEBPACK_IMPORTED_MODULE_1__.OrbitRegime.GEO:
                return [1 + 0.2 * t, 1 + 0.5 * t, 1 + 0.3 * t];
            case _main_js__WEBPACK_IMPORTED_MODULE_1__.OrbitRegime.HEO:
                // Highly elliptical transfer or Molniya
                return [1 + 1.2 * t, 1 + 2.4 * t, 1 + 1.4 * t];
            default:
                return [1 + 0.6 * t, 1 + 1.2 * t, 1 + 0.8 * t];
        }
    }
}


}),
"./src/engine/ootk/src/covariance/StateCovariance.ts": 
/*!***********************************************************!*\
  !*** ./src/engine/ootk/src/covariance/StateCovariance.ts ***!
  \***********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CovarianceFrame: () => (CovarianceFrame),
  StateCovariance: () => (StateCovariance)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/** Covariance Frame */
var CovarianceFrame;
(function (CovarianceFrame) {
    /** Earth-centered inertial */
    CovarianceFrame["ECI"] = "eci";
    /** Radial-Intrack-Crosstrack */
    CovarianceFrame["RIC"] = "ric";
})(CovarianceFrame || (CovarianceFrame = {}));
// / State covariance.
class StateCovariance {
    matrix;
    frame;
    /**
     * Create a new [StateCovariance] object given its covariance [matrix] and
     * [CovarianceFrame].
     * @param matrix The covariance matrix.
     * @param frame The covariance frame.
     * @returns A new [StateCovariance] object.
     */
    constructor(matrix, frame) {
        this.matrix = matrix;
        this.frame = frame;
        // Nothing to do here.
    }
    // / Create a new [StateCovariance] object from 1-sigma values.
    static fromSigmas(sigmas, frame) {
        const n = sigmas.length;
        const output = _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix.zero(n, n);
        for (let i = 0; i < n; i++) {
            output.elements[i][i] = Math.max(sigmas[i] * sigmas[i], 1e-32);
        }
        return new StateCovariance(output, frame);
    }
    /**
     * Calculates the standard deviations (sigmas) of each element in the covariance matrix.
     * @returns A vector containing the standard deviations of each element in the covariance matrix.
     */
    sigmas() {
        const c = this.matrix.columns;
        const result = new Float64Array(c);
        for (let i = 0; i < c; i++) {
            const variance = this.matrix.elements[i][i];
            result[i] = Math.sqrt(variance);
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(result);
    }
}


}),
"./src/engine/ootk/src/covariance/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/covariance/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CovarianceFrame: () => (/* reexport safe */ _StateCovariance_js__WEBPACK_IMPORTED_MODULE_1__.CovarianceFrame),
  CovarianceSample: () => (/* reexport safe */ _CovarianceSample_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceSample),
  StateCovariance: () => (/* reexport safe */ _StateCovariance_js__WEBPACK_IMPORTED_MODULE_1__.StateCovariance)
});
/* ESM import */var _CovarianceSample_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./CovarianceSample.js */ "./src/engine/ootk/src/covariance/CovarianceSample.ts");
/* ESM import */var _StateCovariance_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./StateCovariance.js */ "./src/engine/ootk/src/covariance/StateCovariance.ts");




}),
"./src/engine/ootk/src/data/DataHandler.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/data/DataHandler.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DataHandler: () => (DataHandler)
});
/* ESM import */var _values_Egm96Data_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./values/Egm96Data.js */ "./src/engine/ootk/src/data/values/Egm96Data.ts");
/* ESM import */var _values_HpAtmosphereData_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./values/HpAtmosphereData.js */ "./src/engine/ootk/src/data/values/HpAtmosphereData.ts");
/* ESM import */var _values_Iau1980Data_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./values/Iau1980Data.js */ "./src/engine/ootk/src/data/values/Iau1980Data.ts");
/* ESM import */var _values_LeapSecondData_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./values/LeapSecondData.js */ "./src/engine/ootk/src/data/values/LeapSecondData.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




/**
 * Astrodynamic data management singleton.
 *
 * This class provides access to the various data sets used by the library.
 * It is a singleton class, and can be accessed via the static
 * [getInstance] method.
 */
class DataHandler {
    static instance_ = new DataHandler();
    constructor() {
        // Prevent instantiation.
    }
    /**
     * Returns the singleton instance of the DataHandler class.
     * @returns The singleton instance of the DataHandler class.
     */
    static getInstance() {
        return DataHandler.instance_;
    }
    /**
     * Retrieves the Egm96 coefficients for the given l and m values.
     * @param l - The degree of the coefficient.
     * @param m - The order of the coefficient.
     * @returns The Egm96Entry object containing the coefficients.
     */
    getEgm96Coeffs(l, m) {
        return _values_Egm96Data_js__WEBPACK_IMPORTED_MODULE_0__.egm96Data.getCoeffs(l, m);
    }
    /**
     * Retrieves the IAU 1980 coefficients for the specified row.
     * @param row The row index of the coefficients to retrieve.
     * @returns The IAU 1980 entry containing the coefficients.
     */
    getIau1980Coeffs(row) {
        return _values_Iau1980Data_js__WEBPACK_IMPORTED_MODULE_2__.iau1980Data.getCoeffs(row);
    }
    /**
     * Retrieves the number of leap seconds for a given Julian date.
     * @param jd The Julian date.
     * @returns The number of leap seconds.
     */
    getLeapSeconds(jd) {
        return _values_LeapSecondData_js__WEBPACK_IMPORTED_MODULE_3__.leapSecondData.getLeapSeconds(jd);
    }
    /**
     * Retrieves the atmosphere data for a given height.
     * @param height The height for which to retrieve the atmosphere data.
     * @returns The atmosphere data for the given height, or null if no data is available.
     */
    getHpAtmosphere(height) {
        return _values_HpAtmosphereData_js__WEBPACK_IMPORTED_MODULE_1__.hpAtmosphereData.getAtmosphere(height);
    }
}


}),
"./src/engine/ootk/src/data/values/Egm96Data.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/data/values/Egm96Data.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Egm96Data: () => (Egm96Data),
  egm96Data: () => (egm96Data)
});
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _egm96_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./egm96.js */ "./src/engine/ootk/src/data/values/egm96.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


// / Container for EGM-96 data.
class Egm96Data {
    coeffs_;
    // / Create a new [Egm96Data] object.
    constructor(coeffs) {
        this.coeffs_ = coeffs;
    }
    /**
     * Create a new [Egm96Data] container given a list of EGM-96
     * coefficient tuples [vals].
     * @param vals List of EGM-96 coefficient tuples.
     * @returns A new [Egm96Data] object.
     */
    static fromVals(vals) {
        const output = [];
        for (const v of vals) {
            const [l, m, clm, slm] = v;
            const k = m === 0 ? 1 : 2;
            const a = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_0__.factorial)(l + m);
            const b = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_0__.factorial)(l - m) * (k * (2 * l + 1));
            const nFac = Math.sqrt(a / b);
            const normalizedClm = clm / nFac;
            const normalizedSlm = slm / nFac;
            output.push([l, m, normalizedClm, normalizedSlm]);
        }
        return new Egm96Data(output);
    }
    // / Return de-normalized EGM-96 coefficients for a given [l] and [m] index.
    getCoeffs(l, m) {
        return this.coeffs_[Egm96Data.index_(l, m)];
    }
    // / Return the EGM-96 index for a given [l] and [m] lookup.
    static index_(l, m) {
        return (((l - 2) * (l + 2) + l) >> 1) - 1 + m;
    }
}
// / De-normalized EGM-96 data container.
const egm96Data = Egm96Data.fromVals(_egm96_js__WEBPACK_IMPORTED_MODULE_1__.egm96);


}),
"./src/engine/ootk/src/data/values/HpAtmosphereData.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/data/values/HpAtmosphereData.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  HpAtmosphereData: () => (HpAtmosphereData),
  hpAtmosphereData: () => (hpAtmosphereData)
});
/* ESM import */var _hpAtmosphere_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./hpAtmosphere.js */ "./src/engine/ootk/src/data/values/hpAtmosphere.ts");
/* ESM import */var _HpAtmosphereResult_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./HpAtmosphereResult.js */ "./src/engine/ootk/src/data/values/HpAtmosphereResult.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


// / Container for Harris-Priester atmosphere data.
class HpAtmosphereData {
    table_;
    hMin_;
    hMax_;
    constructor(table) {
        this.table_ = table;
        if (table.length === 0 || typeof table[0]?.[0] === 'undefined') {
            throw new Error('Table must have at least one valid entry.');
        }
        this.hMin_ = table[0][0];
        this.hMax_ = table[table.length - 1]?.[0] ?? 0;
    }
    static fromVals(vals) {
        const output = [];
        for (const v of vals) {
            const [h, minD, maxD] = v;
            output.push([h, minD, maxD]);
        }
        return new HpAtmosphereData(output);
    }
    getAtmosphere(height) {
        if (height < this.hMin_ || height > this.hMax_) {
            return null;
        }
        let index = 0;
        while (index < this.table_.length - 2 && height > (this.table_[index + 1])[0]) {
            index++;
        }
        return new _HpAtmosphereResult_js__WEBPACK_IMPORTED_MODULE_1__.HpAtmosphereResult(height, (this.table_[index]), (this.table_[index + 1]));
    }
}
const hpAtmosphereData = HpAtmosphereData.fromVals(_hpAtmosphere_js__WEBPACK_IMPORTED_MODULE_0__.hpAtmosphere);


}),
"./src/engine/ootk/src/data/values/HpAtmosphereResult.ts": 
/*!***************************************************************!*\
  !*** ./src/engine/ootk/src/data/values/HpAtmosphereResult.ts ***!
  \***************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  HpAtmosphereResult: () => (HpAtmosphereResult)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Harris-Priester atmospheric density bracket.
class HpAtmosphereResult {
    // / Height above Earth's surface _(km)_.
    height;
    // / Lower bound for atmospheric parameters.
    hp0;
    // / Upper bound for atmospheric parameters.
    hp1;
    // / Create a new [HpAtmosphereResult] object.
    constructor(height, hp0, hp1) {
        this.height = height;
        this.hp0 = hp0;
        this.hp1 = hp1;
    }
}


}),
"./src/engine/ootk/src/data/values/Iau1980Data.ts": 
/*!********************************************************!*\
  !*** ./src/engine/ootk/src/data/values/Iau1980Data.ts ***!
  \********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Iau1980Data: () => (Iau1980Data),
  iau1980Data: () => (iau1980Data)
});
/* ESM import */var _iau1980_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./iau1980.js */ "./src/engine/ootk/src/data/values/iau1980.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Container for IAU-1980 data.
class Iau1980Data {
    coeffs_;
    // / Create a new [Iau1980Data] object.
    constructor(coeffs) {
        this.coeffs_ = coeffs;
    }
    /**
     * Create a new [Iau1980Data] container object from an array of IAU-1980
     * coefficient tuples [coeffs].
     * @param coeffs IAU-1980 coefficients.
     * @returns A new [Iau1980Data] object.
     */
    static fromCoeffs(coeffs) {
        const output = [];
        for (const c of coeffs) {
            const [a1, a2, a3, a4, a5, ai, bi, ci, di] = c;
            output.push([a1, a2, a3, a4, a5, ai, bi, ci, di]);
        }
        return new Iau1980Data(output);
    }
    // / Get IAU-1980 coefficients for a given row number.
    getCoeffs(row) {
        return this.coeffs_[row];
    }
}
// / IAU-1980 data container.
const iau1980Data = Iau1980Data.fromCoeffs(_iau1980_js__WEBPACK_IMPORTED_MODULE_0__.iau1980);


}),
"./src/engine/ootk/src/data/values/LeapSecond.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/data/values/LeapSecond.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  LeapSecond: () => (LeapSecond)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Leap second data.
class LeapSecond {
    // / Julian date.
    jd;
    // / Offset seconds.
    offset;
    // / Create a new [LeapSecond] object.
    constructor(jd, offset) {
        this.jd = jd;
        this.offset = offset;
    }
}


}),
"./src/engine/ootk/src/data/values/LeapSecondData.ts": 
/*!***********************************************************!*\
  !*** ./src/engine/ootk/src/data/values/LeapSecondData.ts ***!
  \***********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  leapSecondData: () => (leapSecondData)
});
/* ESM import */var _LeapSecond_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./LeapSecond.js */ "./src/engine/ootk/src/data/values/LeapSecond.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Leap second value tuples.
const leapSeconds = [
    [2441317.5, 10],
    [2441499.5, 11],
    [2441683.5, 12],
    [2442048.5, 13],
    [2442413.5, 14],
    [2442778.5, 15],
    [2443144.5, 16],
    [2443509.5, 17],
    [2443874.5, 18],
    [2444239.5, 19],
    [2444786.5, 20],
    [2445151.5, 21],
    [2445516.5, 22],
    [2446247.5, 23],
    [2447161.5, 24],
    [2447892.5, 25],
    [2448257.5, 26],
    [2448804.5, 27],
    [2449169.5, 28],
    [2449534.5, 29],
    [2450083.5, 30],
    [2450630.5, 31],
    [2451179.5, 32],
    [2453736.5, 33],
    [2454832.5, 34],
    [2456109.5, 35],
    [2457204.5, 36],
    [2457754.5, 37],
];
// / Leap second data container.
class LeapSecondData {
    offsets_;
    jdFirst_;
    jdLast_;
    offsetFirst_;
    offsetLast_;
    constructor(offsets) {
        this.offsets_ = offsets;
        this.jdFirst_ = (this.offsets_[0]).jd;
        this.jdLast_ = (this.offsets_[this.offsets_.length - 1]).jd;
        this.offsetFirst_ = (this.offsets_[0]).offset;
        this.offsetLast_ = (this.offsets_[this.offsets_.length - 1]).offset;
    }
    /**
     * Create a new [LeapSecondData] container given a list of leap second
     * value tuples [vals].
     * @param vals Leap second value tuples.
     * @returns A new [LeapSecondData] container.
     */
    static fromVals(vals) {
        const output = [];
        for (const v of vals) {
            const [jd, offset] = v;
            output.push(new _LeapSecond_js__WEBPACK_IMPORTED_MODULE_0__.LeapSecond(jd, offset));
        }
        return new LeapSecondData(output);
    }
    // / Return the number of leap seconds for a given Julian date [jd].
    getLeapSeconds(jd) {
        if (jd >= this.jdLast_) {
            return this.offsetLast_;
        }
        if (jd <= this.jdFirst_) {
            return this.offsetFirst_;
        }
        for (let i = 0; i < this.offsets_.length - 2; i++) {
            const currentLeapSecond = this.offsets_[i];
            const nextLeapSecond = this.offsets_[i + 1];
            if (jd >= currentLeapSecond.jd && jd < nextLeapSecond.jd) {
                return currentLeapSecond.offset;
            }
        }
        return 0;
    }
}
// / Leap second data container.
const leapSecondData = LeapSecondData.fromVals(leapSeconds);


}),
"./src/engine/ootk/src/data/values/egm96.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/data/values/egm96.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  egm96: () => (egm96)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / The first degree 36 EGM-96 normalized coefficients.
const egm96 = [
    [2, 0, -0.000484165371736, 0],
    [2, 1, -1.86987635955e-10, 1.19528012031e-9],
    [2, 2, 0.00000243914352398, -0.00000140016683654],
    [3, 0, 9.57254173792e-7, 0],
    [3, 1, 0.00000202998882184, 2.48513158716e-7],
    [3, 2, 9.04627768605e-7, -6.19025944205e-7],
    [3, 3, 7.21072657057e-7, 0.00000141435626958],
    [4, 0, 5.39873863789e-7, 0],
    [4, 1, -5.36321616971e-7, -4.73440265853e-7],
    [4, 2, 3.50694105785e-7, 6.6267157254e-7],
    [4, 3, 9.90771803829e-7, -2.00928369177e-7],
    [4, 4, -1.88560802735e-7, 3.08853169333e-7],
    [5, 0, 6.8532347563e-8, 0],
    [5, 1, -6.21012128528e-8, -9.44226127525e-8],
    [5, 2, 6.52438297612e-7, -3.23349612668e-7],
    [5, 3, -4.51955406071e-7, -2.14847190624e-7],
    [5, 4, -2.95301647654e-7, 4.96658876769e-8],
    [5, 5, 1.74971983203e-7, -6.69384278219e-7],
    [6, 0, -1.49957994714e-7, 0],
    [6, 1, -7.60879384947e-8, 2.62890545501e-8],
    [6, 2, 4.81732442832e-8, -3.73728201347e-7],
    [6, 3, 5.71730990516e-8, 9.02694517163e-9],
    [6, 4, -8.62142660109e-8, -4.71408154267e-7],
    [6, 5, -2.6713332549e-7, -5.36488432483e-7],
    [6, 6, 9.67616121092e-9, -2.37192006935e-7],
    [7, 0, 9.0978937145e-8, 0],
    [7, 1, 2.79872910488e-7, 9.54336911867e-8],
    [7, 2, 3.29743816488e-7, 9.30667596042e-8],
    [7, 3, 2.50398657706e-7, -2.17198608738e-7],
    [7, 4, -2.75114355257e-7, -1.23800392323e-7],
    [7, 5, 1.93765507243e-9, 1.77377719872e-8],
    [7, 6, -3.58856860645e-7, 1.51789817739e-7],
    [7, 7, 1.09185148045e-9, 2.44415707993e-8],
    [8, 0, 4.96711667324e-8, 0],
    [8, 1, 2.33422047893e-8, 5.90060493411e-8],
    [8, 2, 8.02978722615e-8, 6.54175425859e-8],
    [8, 3, -1.91877757009e-8, -8.63454445021e-8],
    [8, 4, -2.44600105471e-7, 7.00233016934e-8],
    [8, 5, -2.55352403037e-8, 8.91462164788e-8],
    [8, 6, -6.57361610961e-8, 3.09238461807e-7],
    [8, 7, 6.72811580072e-8, 7.47440473633e-8],
    [8, 8, -1.24092493016e-7, 1.20533165603e-7],
    [9, 0, 2.76714300853e-8, 0],
    [9, 1, 1.43387502749e-7, 2.16834947618e-8],
    [9, 2, 2.22288318564e-8, -3.22196647116e-8],
    [9, 3, -1.60811502143e-7, -7.42287409462e-8],
    [9, 4, -9.00179225336e-9, 1.94666779475e-8],
    [9, 5, -1.66165092924e-8, -5.41113191483e-8],
    [9, 6, 6.26941938248e-8, 2.22903525945e-7],
    [9, 7, -1.18366323475e-7, -9.65152667886e-8],
    [9, 8, 1.88436022794e-7, -3.08566220421e-9],
    [9, 9, -4.77475386132e-8, 9.66412847714e-8],
    [10, 0, 5.26222488569e-8, 0],
    [10, 1, 8.35115775652e-8, -1.31314331796e-7],
    [10, 2, -9.42413882081e-8, -5.1579165739e-8],
    [10, 3, -6.89895048176e-9, -1.53768828694e-7],
    [10, 4, -8.40764549716e-8, -7.92806255331e-8],
    [10, 5, -4.93395938185e-8, -5.05370221897e-8],
    [10, 6, -3.75885236598e-8, -7.95667053872e-8],
    [10, 7, 8.11460540925e-9, -3.36629641314e-9],
    [10, 8, 4.04927981694e-8, -9.18705975922e-8],
    [10, 9, 1.25491334939e-7, -3.76516222392e-8],
    [10, 10, 1.00538634409e-7, -2.4014844952e-8],
    [11, 0, -5.09613707522e-8, 0],
    [11, 1, 1.51687209933e-8, -2.68604146166e-8],
    [11, 2, 1.86309749878e-8, -9.90693862047e-8],
    [11, 3, -3.09871239854e-8, -1.4813180426e-7],
    [11, 4, -3.89580205051e-8, -6.3666651198e-8],
    [11, 5, 3.77848029452e-8, 4.94736238169e-8],
    [11, 6, -1.18676592395e-9, 3.44769584593e-8],
    [11, 7, 4.11565188074e-9, -8.98252808977e-8],
    [11, 8, -5.984108413e-9, 2.43989612237e-8],
    [11, 9, -3.14231072723e-8, 4.17731829829e-8],
    [11, 10, -5.21882681927e-8, -1.83364561788e-8],
    [11, 11, 4.60344448746e-8, -6.96662308185e-8],
    [12, 0, 3.77252636558e-8, 0],
    [12, 1, -5.40654977836e-8, -4.35675748979e-8],
    [12, 2, 1.42979642253e-8, 3.20975937619e-8],
    [12, 3, 3.93995876403e-8, 2.44264863505e-8],
    [12, 4, -6.86908127934e-8, 4.15081109011e-9],
    [12, 5, 3.0941112873e-8, 7.82536279033e-9],
    [12, 6, 3.41523275208e-9, 3.91765484449e-8],
    [12, 7, -1.86909958587e-8, 3.56131849382e-8],
    [12, 8, -2.53769398865e-8, 1.69361024629e-8],
    [12, 9, 4.22880630662e-8, 2.52692598301e-8],
    [12, 10, -6.17619654902e-9, 3.08375794212e-8],
    [12, 11, 1.12502994122e-8, -6.37946501558e-9],
    [12, 12, -2.4953260739e-9, -1.117806019e-8],
    [13, 0, 4.22982206413e-8, 0],
    [13, 1, -5.13569699124e-8, 3.90510386685e-8],
    [13, 2, 5.59217667099e-8, -6.27337565381e-8],
    [13, 3, -2.19360927945e-8, 9.74829362237e-8],
    [13, 4, -3.13762599666e-9, -1.19627874492e-8],
    [13, 5, 5.90049394905e-8, 6.64975958036e-8],
    [13, 6, -3.59038073075e-8, -6.57280613686e-9],
    [13, 7, 2.53002147087e-9, -6.21470822331e-9],
    [13, 8, -9.83150822695e-9, -1.04740222825e-8],
    [13, 9, 2.47325771791e-8, 4.52870369936e-8],
    [13, 10, 4.1032465393e-8, -3.6812102948e-8],
    [13, 11, -4.43869677399e-8, -4.76507804288e-9],
    [13, 12, -3.12622200222e-8, 8.78405809267e-8],
    [13, 13, -6.12759553199e-8, 6.85261488594e-8],
    [14, 0, -2.42786502921e-8, 0],
    [14, 1, -1.86968616381e-8, 2.94747542249e-8],
    [14, 2, -3.67789379502e-8, -5.16779392055e-9],
    [14, 3, 3.58875097333e-8, 2.04618827833e-8],
    [14, 4, 1.83865617792e-9, -2.26780613566e-8],
    [14, 5, 2.87344273542e-8, -1.63882249728e-8],
    [14, 6, -1.94810485574e-8, 2.47831272781e-9],
    [14, 7, 3.75003839415e-8, -4.17291319429e-9],
    [14, 8, -3.50946485865e-8, -1.53515265203e-8],
    [14, 9, 3.20284939341e-8, 2.88804922064e-8],
    [14, 10, 3.90329180008e-8, -1.44308452469e-9],
    [14, 11, 1.53970516502e-8, -3.90548173245e-8],
    [14, 12, 8.40829163869e-9, -3.11327189117e-8],
    [14, 13, 3.22147043964e-8, 4.5189722496e-8],
    [14, 14, -5.18980794309e-8, -4.81506636748e-9],
    [15, 0, 1.47910068708e-9, 0],
    [15, 1, 1.00817268177e-8, 1.09773066324e-8],
    [15, 2, -2.13942673775e-8, -3.08914875777e-8],
    [15, 3, 5.21392929041e-8, 1.72892926103e-8],
    [15, 4, -4.08150084078e-8, 6.50174707794e-9],
    [15, 5, 1.24935723108e-8, 8.08375563996e-9],
    [15, 6, 3.31211643896e-8, -3.68246004304e-8],
    [15, 7, 5.96210699259e-8, 5.31841171879e-9],
    [15, 8, -3.22428691498e-8, 2.21523579587e-8],
    [15, 9, 1.28788268085e-8, 3.75629820829e-8],
    [15, 10, 1.04688722521e-8, 1.47222147015e-8],
    [15, 11, -1.11675061934e-9, 1.80996198432e-8],
    [15, 12, -3.23962134415e-8, 1.55243104746e-8],
    [15, 13, -2.83933019117e-8, -4.22066791103e-9],
    [15, 14, 5.1916885933e-9, -2.43752739666e-8],
    [15, 15, -1.90930538322e-8, -4.71139421558e-9],
    [16, 0, -3.15322986722e-9, 0],
    [16, 1, 2.58360856231e-8, 3.25447560859e-8],
    [16, 2, -2.33671404512e-8, 2.88799363439e-8],
    [16, 3, -3.36019429391e-8, -2.2041898801e-8],
    [16, 4, 4.02316284314e-8, 4.83837716909e-8],
    [16, 5, -1.29501939245e-8, -3.19458578129e-9],
    [16, 6, 1.40239252323e-8, -3.50760208303e-8],
    [16, 7, -7.08412635136e-9, -8.81581561131e-9],
    [16, 8, -2.09018868094e-8, 5.0052739053e-9],
    [16, 9, -2.18588720643e-8, -3.95012419994e-8],
    [16, 10, -1.17529900814e-8, 1.14211582961e-8],
    [16, 11, 1.87574042592e-8, -3.03161919925e-9],
    [16, 12, 1.95400194038e-8, 6.66983574071e-9],
    [16, 13, 1.38196369576e-8, 1.02778499508e-9],
    [16, 14, -1.93182168856e-8, -3.86174893776e-8],
    [16, 15, -1.45149060142e-8, -3.27443078739e-8],
    [16, 16, -3.79671710746e-8, 3.02155372655e-9],
    [17, 0, 1.97605066395e-8, 0],
    [17, 1, -2.54177575118e-8, -3.06630529689e-8],
    [17, 2, -1.95988656721e-8, 6.4926589341e-9],
    [17, 3, 5.64123066224e-9, 6.78327095529e-9],
    [17, 4, 7.07457075637e-9, 2.49437600834e-8],
    [17, 5, -1.54987006052e-8, 6.60021551851e-9],
    [17, 6, -1.18194012847e-8, -2.89770975177e-8],
    [17, 7, 2.42149702381e-8, -4.22222973697e-9],
    [17, 8, 3.88442097559e-8, 3.58904095943e-9],
    [17, 9, 3.81356493231e-9, -2.81466943714e-8],
    [17, 10, -3.88216085542e-9, 1.81328176508e-8],
    [17, 11, -1.57356600363e-8, 1.06560649404e-8],
    [17, 12, 2.88013010655e-8, 2.03450136084e-8],
    [17, 13, 1.65503425731e-8, 2.04667531435e-8],
    [17, 14, -1.41983872649e-8, 1.14948025244e-8],
    [17, 15, 5.42100361657e-9, 5.32610369811e-9],
    [17, 16, -3.01992205043e-8, 3.65331918531e-9],
    [17, 17, -3.43086856041e-8, -1.98523455381e-8],
    [18, 0, 5.08691038332e-9, 0],
    [18, 1, 7.21098449649e-9, -3.88714473013e-8],
    [18, 2, 1.40631771205e-8, 1.00093396253e-8],
    [18, 3, -5.07232520873e-9, -4.90865931335e-9],
    [18, 4, 5.48759308217e-8, -1.3526711772e-9],
    [18, 5, 5.48710485555e-9, 2.64338629459e-8],
    [18, 6, 1.46570755271e-8, -1.36438019951e-8],
    [18, 7, 6.75812328417e-9, 6.88577494235e-9],
    [18, 8, 3.07619845144e-8, 4.17827734107e-9],
    [18, 9, -1.8847060188e-8, 3.68302736953e-8],
    [18, 10, 5.27535358934e-9, -4.66091535881e-9],
    [18, 11, -7.2962851896e-9, 1.9521520802e-9],
    [18, 12, -2.97449412422e-8, -1.64497878395e-8],
    [18, 13, -6.27919717152e-9, -3.48383939938e-8],
    [18, 14, -8.1560533641e-9, -1.28636585027e-8],
    [18, 15, -4.05003412879e-8, -2.02684998021e-8],
    [18, 16, 1.04141042028e-8, 6.61468817624e-9],
    [18, 17, 3.58771586841e-9, 4.48065587564e-9],
    [18, 18, 3.12351953717e-9, -1.09906032543e-8],
    [19, 0, -3.25780965394e-9, 0],
    [19, 1, -7.59903885319e-9, 1.26835472605e-9],
    [19, 2, 3.53541528655e-8, -1.31346303514e-9],
    [19, 3, -9.74103607309e-9, 1.50662259043e-9],
    [19, 4, 1.57039009057e-8, -7.61677383811e-9],
    [19, 5, 1.09629213379e-8, 2.83172176438e-8],
    [19, 6, -4.08745178658e-9, 1.86219430719e-8],
    [19, 7, 4.78275337044e-9, -7.172834559e-9],
    [19, 8, 2.9490836428e-8, -9.93037002883e-9],
    [19, 9, 3.07961427159e-9, 6.94110477214e-9],
    [19, 10, -3.38415069043e-8, -7.37981767136e-9],
    [19, 11, 1.60443652916e-8, 9.96673453483e-9],
    [19, 12, -2.47106581581e-9, 9.16852310642e-9],
    [19, 13, -7.4471737998e-9, -2.82584466742e-8],
    [19, 14, -4.70502589215e-9, -1.29526697983e-8],
    [19, 15, -1.76580549771e-8, -1.40350990039e-8],
    [19, 16, -2.16950096188e-8, -7.24534721567e-9],
    [19, 17, 2.90444936079e-8, -1.5345653107e-8],
    [19, 18, 3.48382199593e-8, -9.54146344917e-9],
    [19, 19, -2.5734934943e-9, 4.83151822363e-9],
    [20, 0, 2.22384610651e-8, 0],
    [20, 1, 5.16303125218e-9, 6.69626726966e-9],
    [20, 2, 1.98831128238e-8, 1.75183843257e-8],
    [20, 3, -3.62601436785e-9, 3.79590724141e-8],
    [20, 4, 2.42238118652e-9, -2.11057611874e-8],
    [20, 5, -1.07042562564e-8, -7.71860083169e-9],
    [20, 6, 1.1047483757e-8, -2.17720365898e-9],
    [20, 7, -2.10090282728e-8, -2.23491503969e-11],
    [20, 8, 4.42419185637e-9, 1.83035804593e-9],
    [20, 9, 1.78846216942e-8, -6.63940865358e-9],
    [20, 10, -3.25394919988e-8, -5.12308873621e-9],
    [20, 11, 1.38992707697e-8, -1.87706454942e-8],
    [20, 12, -6.3575060075e-9, 1.80260853103e-8],
    [20, 13, 2.75222725997e-8, 6.90887077588e-9],
    [20, 14, 1.15841169405e-8, -1.43176160143e-8],
    [20, 15, -2.60130744291e-8, -7.84379672413e-10],
    [20, 16, -1.24137147118e-8, -2.77500443628e-10],
    [20, 17, 4.3690966796e-9, -1.37420446198e-8],
    [20, 18, 1.51842883022e-8, -8.08429903142e-10],
    [20, 19, -3.14942002852e-9, 1.06505202245e-8],
    [20, 20, 4.01448327968e-9, -1.20450644785e-8],
    [21, 0, 5.87820252575e-9, 0],
    [21, 1, -1.61000670141e-8, 2.84359400791e-8],
    [21, 2, -6.54460482558e-9, 3.78474868508e-9],
    [21, 3, 1.9549199526e-8, 2.26286963716e-8],
    [21, 4, -5.76604339239e-9, 1.94493782631e-8],
    [21, 5, 2.58856303016e-9, 1.70850368669e-9],
    [21, 6, -1.40168810589e-8, -2.73814826381e-12],
    [21, 7, -8.64357168475e-9, 4.42612277119e-9],
    [21, 8, -1.70477278237e-8, 1.5071119263e-9],
    [21, 9, 1.64489062394e-8, 8.30113196365e-9],
    [21, 10, -1.09928976409e-8, -1.46913794684e-9],
    [21, 11, 6.99300364214e-9, -3.53590565124e-8],
    [21, 12, -3.19300109594e-9, 1.45786917947e-8],
    [21, 13, -1.8985452459e-8, 1.40514791436e-8],
    [21, 14, 2.03580785674e-8, 7.5577246284e-9],
    [21, 15, 1.75530220278e-8, 1.04533886832e-8],
    [21, 16, 7.86969109367e-9, -6.56089715279e-9],
    [21, 17, -6.99484489981e-9, -7.36064901147e-9],
    [21, 18, 2.59643291521e-8, -1.1156080613e-8],
    [21, 19, -2.7374163641e-8, 1.63958190052e-8],
    [21, 20, -2.68682473584e-8, 1.62086057168e-8],
    [21, 21, 8.30374873932e-9, -3.75546121742e-9],
    [22, 0, -1.13735124259e-8, 0],
    [22, 1, 1.62309865679e-8, -3.77303475153e-9],
    [22, 2, -2.64090261387e-8, -2.10832402428e-9],
    [22, 3, 1.1658001654e-8, 1.06764617222e-8],
    [22, 4, -2.70979141451e-9, 1.74980820565e-8],
    [22, 5, -1.8645262501e-9, 7.44718166476e-10],
    [22, 6, 9.64390704406e-9, -6.37316743908e-9],
    [22, 7, 1.59715981795e-8, 4.39600942993e-9],
    [22, 8, -2.35157426998e-8, 4.83673695086e-9],
    [22, 9, 8.29435796737e-9, 8.73382159986e-9],
    [22, 10, 6.00704037701e-9, 2.21854121109e-8],
    [22, 11, -4.96078301539e-9, -1.78822672474e-8],
    [22, 12, 2.13502315463e-9, -7.96120522503e-9],
    [22, 13, -1.72631843979e-8, 1.97026896892e-8],
    [22, 14, 1.09297133018e-8, 8.25280905301e-9],
    [22, 15, 2.58410840629e-8, 4.60172998318e-9],
    [22, 16, 1.41258558921e-10, -7.182380053e-9],
    [22, 17, 8.89294096846e-9, -1.45618348246e-8],
    [22, 18, 1.05047447464e-8, -1.64271275481e-8],
    [22, 19, 1.41305509124e-8, -3.84537168599e-9],
    [22, 20, -1.67617655441e-8, 1.99561513321e-8],
    [22, 21, -2.50948756455e-8, 2.36151346133e-8],
    [22, 22, -9.59596694809e-9, 2.49861413883e-9],
    [23, 0, -2.26201075082e-8, 0],
    [23, 1, 1.10870239758e-8, 1.6137915153e-8],
    [23, 2, -1.35191027779e-8, -5.01411714852e-9],
    [23, 3, -2.45128011445e-8, -1.60570438998e-8],
    [23, 4, -2.39887874558e-8, 7.31536362289e-9],
    [23, 5, 7.99636624146e-10, -1.6144974141e-10],
    [23, 6, -1.26082781309e-8, 1.61308155632e-8],
    [23, 7, -8.04132133762e-9, -1.11647197494e-9],
    [23, 8, 7.53785326469e-9, -3.2967992522e-10],
    [23, 9, 2.5505325495e-9, -1.28071525548e-8],
    [23, 10, 1.65167929134e-8, -1.85239620853e-9],
    [23, 11, 9.42656822725e-9, 1.52386181583e-8],
    [23, 12, 1.63632625535e-8, -1.24098327824e-8],
    [23, 13, -1.15107832808e-8, -4.84279171627e-9],
    [23, 14, 6.75321602206e-9, -1.82899962212e-9],
    [23, 15, 1.8689804286e-8, -3.60523754481e-9],
    [23, 16, 6.13840121864e-9, 1.10362707266e-8],
    [23, 17, -5.5372102391e-9, -1.2845906046e-8],
    [23, 18, 8.43361263813e-9, -1.49115921605e-8],
    [23, 19, -5.20848228342e-9, 1.07789593943e-8],
    [23, 20, 8.60434396837e-9, -5.34641639372e-9],
    [23, 21, 1.54578189867e-8, 1.15333325358e-8],
    [23, 22, -1.78417206471e-8, 4.33092348903e-9],
    [23, 23, 2.85393980111e-9, -1.1323294597e-8],
    [24, 0, 7.63657386411e-10, 0],
    [24, 1, -3.14943681427e-9, -1.77191190396e-9],
    [24, 2, 1.38595572093e-9, 1.711040664e-8],
    [24, 3, -4.76406913528e-9, -9.42329378125e-9],
    [24, 4, 6.05108036341e-9, 5.49769910191e-9],
    [24, 5, -7.2947904748e-9, -2.13826490504e-8],
    [24, 6, 4.54210367535e-9, 1.85596665318e-9],
    [24, 7, -6.14244489298e-9, 4.70081667951e-9],
    [24, 8, 1.54822444425e-8, -4.34472097787e-9],
    [24, 9, -9.76623425797e-9, -1.6275513762e-8],
    [24, 10, 1.08934628974e-8, 2.09168783608e-8],
    [24, 11, 1.45280775337e-8, 1.87398018797e-8],
    [24, 12, 1.18970310717e-8, -6.2293309815e-9],
    [24, 13, -2.89676673058e-9, 3.13251295024e-9],
    [24, 14, -2.00006558603e-8, -1.87249636821e-9],
    [24, 15, 6.10396350698e-9, -1.58957680563e-8],
    [24, 16, 8.88750753375e-9, 2.96492703352e-9],
    [24, 17, -1.19629964611e-8, -5.82074593955e-9],
    [24, 18, -6.52630641555e-10, -1.01332355837e-8],
    [24, 19, -4.38896550264e-9, -8.14552569977e-9],
    [24, 20, -5.17551981851e-9, 8.90354942378e-9],
    [24, 21, 6.03436755046e-9, 1.40116090741e-8],
    [24, 22, 3.93640283055e-9, -4.28327655754e-9],
    [24, 23, -6.1428347955e-9, -8.692679021e-9],
    [24, 24, 1.23903921309e-8, -3.75059286959e-9],
    [25, 0, 3.21309208115e-9, 0],
    [25, 1, 6.89649208567e-9, -7.995518294e-9],
    [25, 2, 2.19498139173e-8, 9.01370249111e-9],
    [25, 3, -1.17774931587e-8, -1.26719024392e-8],
    [25, 4, 9.4254362892e-9, 6.84937199311e-10],
    [25, 5, -1.00497487339e-8, -9.2212239967e-10],
    [25, 6, 1.66832871654e-8, 4.30583576199e-10],
    [25, 7, 7.71426681671e-9, -4.11703290425e-9],
    [25, 8, 3.1565194415e-9, -7.81960217669e-10],
    [25, 9, -2.99385350515e-8, 2.12695473199e-8],
    [25, 10, 8.81931818034e-9, -4.18041586166e-9],
    [25, 11, 1.2340148568e-9, 1.08069128123e-8],
    [25, 12, -7.65146786755e-9, 1.1747374286e-8],
    [25, 13, 8.32308127158e-9, -1.13072604626e-8],
    [25, 14, -1.97042124794e-8, 6.53183488635e-9],
    [25, 15, -4.35732052985e-9, -7.35147227573e-9],
    [25, 16, 9.18239548455e-10, -1.28124888592e-8],
    [25, 17, -1.52176535379e-8, -3.21280397924e-9],
    [25, 18, 1.21901534245e-9, -1.49040483259e-8],
    [25, 19, 7.77589111757e-9, 9.92518771941e-9],
    [25, 20, -7.50856670672e-9, -5.62826155305e-10],
    [25, 21, 1.0723284068e-8, 8.16090174381e-9],
    [25, 22, -1.39902235929e-8, 3.58546198324e-9],
    [25, 23, 8.40270853655e-9, -1.23338407961e-8],
    [25, 24, 4.12447134569e-9, -8.30716465317e-9],
    [25, 25, 1.07484366767e-8, 4.72369913984e-9],
    [26, 0, 5.05833635414e-9, 0],
    [26, 1, -1.54756177965e-9, -7.70012788871e-9],
    [26, 2, -3.58729876836e-9, 1.14484111182e-8],
    [26, 3, 1.40505671267e-8, 4.30905534294e-9],
    [26, 4, 1.90548709216e-8, -1.94161179658e-8],
    [26, 5, 1.07190025408e-8, 9.08952851813e-9],
    [26, 6, 1.13116909406e-8, -9.34393384449e-9],
    [26, 7, -1.562282956e-9, 4.81168302477e-9],
    [26, 8, 3.94920146317e-9, 1.153405253e-9],
    [26, 9, -1.20371433638e-8, 4.75177058134e-10],
    [26, 10, -1.41246124334e-8, -6.45217247294e-9],
    [26, 11, -5.20385857649e-9, 2.12443340407e-9],
    [26, 12, -1.75071176484e-8, 2.01974971938e-9],
    [26, 13, -3.35708835245e-11, 1.50474091686e-9],
    [26, 14, 7.96385051492e-9, 7.84704068835e-9],
    [26, 15, -1.32388781089e-8, 8.03960091442e-9],
    [26, 16, 1.29093226253e-9, -6.11434455706e-9],
    [26, 17, -1.24494157564e-8, 7.8077484564e-9],
    [26, 18, -1.30317424459e-8, 4.9998916257e-9],
    [26, 19, -2.05807464595e-9, 3.54396135438e-9],
    [26, 20, 6.55952144018e-9, -1.1687804118e-8],
    [26, 21, -8.70038868454e-9, 1.68222257564e-9],
    [26, 22, 1.01580452049e-8, 7.54358531576e-9],
    [26, 23, 1.24105057436e-9, 1.08580088935e-8],
    [26, 24, 8.58620351967e-9, 1.48288510099e-8],
    [26, 25, 3.93441578873e-9, -5.97792415806e-10],
    [26, 26, 3.93179749568e-10, 1.93894997772e-9],
    [27, 0, 2.7717632236e-9, 0],
    [27, 1, 2.48982909452e-9, 3.77378455357e-9],
    [27, 2, 1.45270146453e-9, 5.03113268026e-10],
    [27, 3, -3.62306812856e-10, 1.088457625e-8],
    [27, 4, -5.99191537157e-10, 9.40517681233e-9],
    [27, 5, 1.67690560888e-8, 1.38338587209e-8],
    [27, 6, 3.64265989803e-9, 6.13032807744e-9],
    [27, 7, -1.23459266009e-8, -3.86514075952e-9],
    [27, 8, -6.1040764482e-9, -8.99504471581e-9],
    [27, 9, 3.40113157078e-9, 1.10992938665e-8],
    [27, 10, -1.33158893187e-8, 1.72832279915e-10],
    [27, 11, 1.98322808107e-9, -9.69054254426e-9],
    [27, 12, -1.13695413044e-8, 1.90072943781e-9],
    [27, 13, -4.97224781272e-9, -4.14521559996e-9],
    [27, 14, 1.55033957088e-8, 1.1882128969e-8],
    [27, 15, -1.80057326196e-9, 1.1763698622e-9],
    [27, 16, 2.7572995289e-9, 2.78770269194e-9],
    [27, 17, 3.79281571763e-9, 3.14983101049e-10],
    [27, 18, -2.87144071715e-9, 7.44190558718e-9],
    [27, 19, -3.26518614707e-10, -2.93243500455e-9],
    [27, 20, -8.55182561846e-10, 3.47617208115e-9],
    [27, 21, 4.86877030983e-9, -7.0872528354e-9],
    [27, 22, -5.74332100084e-9, 2.90056687384e-9],
    [27, 23, -5.41033470941e-9, -1.10452433655e-8],
    [27, 24, 4.16951885933e-10, -1.80038186307e-9],
    [27, 25, 1.22815470212e-8, 5.62425137285e-9],
    [27, 26, -6.59498075164e-9, -2.22838418639e-9],
    [27, 27, 7.60067381059e-9, 6.9238741892e-10],
    [28, 0, -9.10376375863e-9, 0],
    [28, 1, -5.55484993587e-9, 7.9330019258e-9],
    [28, 2, -1.5189131211e-8, -7.97957089012e-9],
    [28, 3, 2.5318254224e-9, 1.11373049392e-8],
    [28, 4, -1.99212752126e-9, 1.25054704704e-8],
    [28, 5, 1.08871875702e-8, -4.22573826989e-9],
    [28, 6, -5.22194316032e-9, 1.32656509709e-8],
    [28, 7, -7.05588863746e-10, 5.12740997711e-9],
    [28, 8, -4.23704976329e-9, -3.32584474553e-9],
    [28, 9, 1.13842461859e-8, -1.04163010811e-8],
    [28, 10, -9.22867885082e-9, 8.17851851593e-9],
    [28, 11, -2.9809734257e-9, -1.45944538949e-9],
    [28, 12, -4.83471863256e-10, 9.64951845027e-9],
    [28, 13, 1.64993974957e-9, 6.63803768689e-9],
    [28, 14, -8.23334828619e-9, -1.26939492243e-8],
    [28, 15, -1.22774798187e-8, -1.97537366262e-9],
    [28, 16, -3.57280690709e-9, -1.35890044766e-8],
    [28, 17, 1.33742628184e-8, -4.72374226319e-9],
    [28, 18, 5.62532322748e-9, -3.87230727328e-9],
    [28, 19, 5.77104709635e-9, 2.35011734292e-8],
    [28, 20, -1.15922189521e-9, 6.62939940662e-9],
    [28, 21, 6.63154344375e-9, 6.33201211223e-9],
    [28, 22, -1.94231451662e-9, -7.33725263107e-9],
    [28, 23, 6.20158165102e-9, 2.61202437682e-9],
    [28, 24, 1.11186270621e-8, -1.35606378769e-8],
    [28, 25, 7.29495896149e-9, -1.76041477031e-8],
    [28, 26, 1.23084992259e-8, 3.89251843939e-9],
    [28, 27, -8.11971206724e-9, 1.3027922855e-9],
    [28, 28, 6.9872587832e-9, 6.80526167979e-9],
    [29, 0, -4.97406439473e-9, 0],
    [29, 1, 4.98979084585e-9, -9.82512461189e-9],
    [29, 2, -3.12119754621e-9, -2.63433487676e-9],
    [29, 3, 1.82518120454e-9, -1.05769977751e-8],
    [29, 4, -2.42786314995e-8, 2.26110758622e-9],
    [29, 5, -6.8110306367e-9, 6.01242555817e-9],
    [29, 6, 1.19592879211e-8, 9.7020069574e-9],
    [29, 7, -5.91100934209e-9, -2.14599788734e-9],
    [29, 8, -1.6946723555e-8, 1.11160276839e-8],
    [29, 9, -1.2937116169e-9, 1.41793573226e-9],
    [29, 10, 1.37184624798e-8, 1.79543486167e-9],
    [29, 11, -5.96272885876e-9, 6.33350180946e-9],
    [29, 12, -4.56278910357e-10, -5.01222008898e-9],
    [29, 13, -1.09095923049e-9, -2.34179014389e-9],
    [29, 14, -3.23718965114e-9, -4.58306325034e-9],
    [29, 15, -9.57359749406e-9, -6.77546725808e-9],
    [29, 16, 1.37450063496e-9, -1.4864526654e-8],
    [29, 17, -1.57662415501e-9, -3.92506699434e-9],
    [29, 18, -3.67597840865e-9, -2.58549575294e-9],
    [29, 19, -6.30046143533e-9, 5.86840708296e-9],
    [29, 20, -7.96446331531e-9, 5.74239983127e-9],
    [29, 21, -9.8726430286e-9, -5.51700601596e-9],
    [29, 22, 1.15574836058e-8, -1.47663300854e-9],
    [29, 23, -1.84576717899e-9, 2.63546763516e-9],
    [29, 24, 3.42199668119e-10, -2.38230581193e-9],
    [29, 25, 5.85864038329e-9, 8.68333958543e-9],
    [29, 26, 7.87039835357e-9, -6.92232980921e-9],
    [29, 27, -7.98313300841e-9, -1.01903214091e-9],
    [29, 28, 9.73355537526e-9, -5.71293958601e-9],
    [29, 29, 1.28224843767e-8, -5.01548480482e-9],
    [30, 0, 6.02882084759e-9, 0],
    [30, 1, -5.57556615596e-10, 1.24285275602e-9],
    [30, 2, -1.0370644769e-8, -2.61802322444e-9],
    [30, 3, 2.14692300603e-9, -1.36464188501e-8],
    [30, 4, -4.55090433473e-10, -3.91117213505e-9],
    [30, 5, -4.36973977446e-9, -5.35558974983e-9],
    [30, 6, 3.28451285815e-10, 3.17808233981e-9],
    [30, 7, 4.04923220309e-9, 1.83962458779e-9],
    [30, 8, 2.54952865236e-9, 4.62058281854e-9],
    [30, 9, -7.32592511128e-9, -9.7277817424e-9],
    [30, 10, 4.27609484555e-9, -4.10864961814e-9],
    [30, 11, -1.04043005227e-8, 1.07581457651e-8],
    [30, 12, 1.71622295302e-8, -1.08456775556e-8],
    [30, 13, 1.42173587056e-8, 2.96806226352e-9],
    [30, 14, 5.11505860834e-9, 8.07288811257e-9],
    [30, 15, 2.10512146846e-10, -1.04541123836e-9],
    [30, 16, -1.08921920457e-8, 4.35254063533e-9],
    [30, 17, -6.14382436271e-9, -6.03140938575e-9],
    [30, 18, -1.1114926509e-8, -7.65521957976e-9],
    [30, 19, -1.2967398433e-8, 2.42005669694e-9],
    [30, 20, -4.89261172033e-9, 1.27655684422e-8],
    [30, 21, -1.0628473781e-8, -5.97537587412e-9],
    [30, 22, -4.83763240001e-9, -9.37720111156e-9],
    [30, 23, 5.7411388543e-9, -1.03756082222e-8],
    [30, 24, -2.35238020789e-9, -2.7590933962e-9],
    [30, 25, 3.04426404856e-9, -1.54853389229e-8],
    [30, 26, 1.22149787623e-9, 1.24069551653e-8],
    [30, 27, -7.95063844863e-9, 1.27529431593e-8],
    [30, 28, -5.47120800289e-9, -7.96006293513e-9],
    [30, 29, 4.1592295424e-9, 1.89489104417e-9],
    [30, 30, 2.64794018006e-9, 8.12994755178e-9],
    [31, 0, 7.33100089318e-9, 0],
    [31, 1, 6.11169376734e-9, -1.60774540844e-8],
    [31, 2, 7.49625106123e-9, 6.37776322444e-9],
    [31, 3, -8.89920966189e-9, -7.6550294416e-9],
    [31, 4, 1.22555580723e-8, -4.94466436575e-9],
    [31, 5, -8.71279064045e-9, 3.08325747379e-9],
    [31, 6, -1.68890803585e-9, 1.3703621527e-9],
    [31, 7, -2.71996133536e-9, -6.8862512168e-10],
    [31, 8, -7.50260355354e-10, 2.28102724239e-9],
    [31, 9, -6.55840403272e-10, 5.24179002617e-9],
    [31, 10, 3.99161675027e-9, -4.73500202132e-9],
    [31, 11, 6.93506892777e-10, 2.08668068881e-8],
    [31, 12, 5.5287540984e-10, 4.52042167068e-9],
    [31, 13, 9.40389423562e-9, 4.6684078573e-9],
    [31, 14, -7.88650771167e-9, 3.51952460147e-9],
    [31, 15, 4.29954776132e-9, -2.80870684394e-9],
    [31, 16, -7.19430261173e-9, 6.11805049979e-9],
    [31, 17, -2.53821168958e-9, 6.83008216722e-9],
    [31, 18, -6.02099321996e-10, -2.04187286905e-9],
    [31, 19, 2.89086482301e-9, 4.43976791609e-9],
    [31, 20, -1.75732193914e-9, 5.64081954558e-9],
    [31, 21, -9.67143669208e-9, 7.09357408027e-9],
    [31, 22, -9.0531201252e-9, -1.18308417466e-8],
    [31, 23, 8.32234353898e-9, 4.51774572555e-9],
    [31, 24, -2.81565064366e-9, -3.34369513768e-9],
    [31, 25, -1.64574268169e-8, -2.20460908971e-9],
    [31, 26, -1.26653070356e-8, 1.59189398991e-9],
    [31, 27, -1.34953305827e-9, 1.07507650019e-8],
    [31, 28, 1.04226918411e-8, 2.8072229491e-9],
    [31, 29, -1.5812688103e-9, -2.18247510672e-9],
    [31, 30, -9.47416722001e-10, -7.78077525656e-9],
    [31, 31, -8.59193452715e-9, -1.85200316483e-9],
    [32, 0, -2.33966288032e-9, 0],
    [32, 1, -1.69210486076e-9, 1.27760467976e-9],
    [32, 2, 1.13999662663e-8, -3.35609127916e-9],
    [32, 3, -1.444433154e-10, 4.05424830941e-9],
    [32, 4, 8.56367829112e-10, -6.75422476107e-9],
    [32, 5, 8.60776205333e-9, 1.82572279646e-9],
    [32, 6, -1.00402568672e-8, -7.6305617634e-9],
    [32, 7, 1.37058613278e-9, 2.75465347035e-9],
    [32, 8, 1.19653531908e-8, 4.91018212548e-9],
    [32, 9, 7.332252213e-9, 7.18971591052e-10],
    [32, 10, 9.12133506379e-11, -5.70680927495e-9],
    [32, 11, -5.42043742127e-9, 7.583606425e-9],
    [32, 12, -1.70289059214e-8, 1.40808168623e-8],
    [32, 13, 4.02186822027e-9, 5.34936491964e-9],
    [32, 14, -5.44420334437e-9, 2.20410694316e-9],
    [32, 15, 5.1658020828e-9, -8.74727531741e-9],
    [32, 16, 4.14867061294e-9, 4.27270420004e-9],
    [32, 17, -6.46857778906e-9, 1.01916486215e-8],
    [32, 18, 1.27286345117e-8, -1.12136888089e-9],
    [32, 19, 7.55189536923e-10, -2.7754653073e-9],
    [32, 20, 3.8161056442e-9, 3.19534855653e-10],
    [32, 21, -2.33262996771e-9, 1.16411650251e-8],
    [32, 22, -1.20880678762e-8, -2.72691793232e-9],
    [32, 23, 8.18682122143e-9, -2.33549712722e-9],
    [32, 24, -3.55036315667e-9, 6.54834763861e-10],
    [32, 25, -1.89374992503e-8, -6.43429532848e-9],
    [32, 26, 5.22535531492e-9, -3.68856221241e-9],
    [32, 27, -4.53740085214e-9, -6.68075560111e-9],
    [32, 28, 1.653041745e-9, -5.73130340772e-9],
    [32, 29, 4.32768192965e-9, 2.88179889934e-9],
    [32, 30, -6.74805866294e-9, 1.39346268546e-9],
    [32, 31, -6.26740251766e-9, -2.18475608171e-10],
    [32, 32, 3.3975660331e-9, 1.42646165155e-9],
    [33, 0, -3.49357179498e-9, 0],
    [33, 1, -1.39642913445e-9, -2.16391760811e-9],
    [33, 2, -7.48774194896e-9, -5.0187208152e-10],
    [33, 3, -1.99661955793e-9, 7.0930410268e-9],
    [33, 4, -4.270199819e-9, 2.27426656698e-9],
    [33, 5, 2.37784729729e-10, 3.74439169451e-9],
    [33, 6, 1.22603039921e-9, -2.87328300836e-9],
    [33, 7, -6.11215086076e-9, 2.49383366316e-9],
    [33, 8, -8.23144405057e-10, 1.44915555407e-8],
    [33, 9, 5.05097392033e-9, 7.4051746902e-9],
    [33, 10, -2.39709923317e-9, 1.07022906758e-9],
    [33, 11, 2.43388836443e-9, -8.67071813487e-9],
    [33, 12, -2.33510532329e-9, 8.9435069891e-9],
    [33, 13, 2.6041538193e-9, 3.13805750981e-9],
    [33, 14, 4.92959662302e-9, 5.71204550617e-9],
    [33, 15, -4.64145303396e-9, -3.47835302325e-9],
    [33, 16, 7.39530517571e-9, 6.28613189283e-9],
    [33, 17, -5.73064590551e-9, 1.28779114927e-8],
    [33, 18, -9.74285933562e-9, -1.89598124592e-9],
    [33, 19, 8.52447331156e-9, 2.07561717246e-9],
    [33, 20, -3.32627500309e-9, -7.77689999053e-9],
    [33, 21, 9.38761672387e-10, 8.17787598674e-10],
    [33, 22, -1.05439940875e-8, -1.56190227392e-8],
    [33, 23, 1.15896250314e-10, -1.01356350767e-8],
    [33, 24, 1.11416074527e-8, -8.57153776484e-9],
    [33, 25, 5.24730532375e-9, -1.04941656537e-8],
    [33, 26, 1.09590005596e-8, 4.5404144025e-9],
    [33, 27, -1.32772908147e-9, 1.26154161942e-9],
    [33, 28, 1.75943381421e-9, -1.02060346415e-9],
    [33, 29, -1.63075128633e-8, 5.72191328891e-9],
    [33, 30, -1.56977064277e-9, -1.84579402264e-8],
    [33, 31, 4.69481868853e-9, 1.02290050028e-9],
    [33, 32, 6.56775919022e-9, -4.39711913398e-9],
    [33, 33, -1.52043850303e-9, 8.31263004529e-9],
    [34, 0, -9.08833340447e-9, 0],
    [34, 1, -2.76889795047e-9, 6.3891897021e-9],
    [34, 2, 6.7688190654e-9, 5.30082118696e-9],
    [34, 3, 1.25429669786e-8, 8.11619669834e-9],
    [34, 4, -8.30005417504e-9, 1.19586870272e-9],
    [34, 5, -3.88131685638e-9, 3.54963449977e-9],
    [34, 6, 4.84093709579e-10, 7.62975480293e-9],
    [34, 7, 2.75125793239e-9, -6.56263573163e-9],
    [34, 8, -9.83446807592e-9, 4.68751478021e-9],
    [34, 9, 1.53042494664e-9, 2.10165697829e-9],
    [34, 10, -7.52633242389e-9, 1.46544229781e-9],
    [34, 11, -3.82043431506e-9, -1.07829735599e-9],
    [34, 12, 1.42629362262e-8, -4.60063642968e-9],
    [34, 13, -3.56240984255e-9, 1.03329523096e-9],
    [34, 14, -2.50187664392e-9, 9.64686908241e-9],
    [34, 15, 3.75939804157e-10, 6.2628624977e-9],
    [34, 16, -1.45874042713e-9, -1.4938092908e-9],
    [34, 17, -4.73747570512e-9, 3.93698829389e-9],
    [34, 18, -1.47488701345e-8, -5.38197998817e-9],
    [34, 19, -3.59837568897e-9, 7.15302015583e-9],
    [34, 20, 3.64466859655e-9, -1.01824147346e-8],
    [34, 21, -9.81980297066e-10, -7.42166456548e-9],
    [34, 22, -3.18152215406e-9, 3.36620175035e-9],
    [34, 23, -1.1297312057e-9, -1.18981902172e-8],
    [34, 24, 8.78079044954e-9, 4.20436158037e-9],
    [34, 25, 8.41097170248e-9, -9.86300815266e-9],
    [34, 26, 3.99964384231e-9, -1.29360014691e-8],
    [34, 27, 1.31566196208e-8, -3.91137836409e-9],
    [34, 28, -1.65320604713e-10, -2.00370653858e-8],
    [34, 29, 7.08151676681e-9, -4.31563574113e-9],
    [34, 30, -2.05666035677e-8, -5.86948946952e-10],
    [34, 31, -4.57411268111e-9, -1.60852780125e-9],
    [34, 32, 9.14033593474e-9, 2.31645138264e-9],
    [34, 33, 1.37617937967e-8, 4.3547198646e-9],
    [34, 34, -8.54011998155e-9, 1.65364599023e-9],
    [35, 0, 8.60443158492e-9, 0],
    [35, 1, -1.07631176168e-8, -1.03576288219e-8],
    [35, 2, -1.48166749807e-8, 7.47316845223e-9],
    [35, 3, 1.88623900305e-9, 3.49967679465e-9],
    [35, 4, -2.82338523108e-9, 9.20674937921e-9],
    [35, 5, -7.23688443416e-9, -1.15478796146e-8],
    [35, 6, 3.28708320436e-9, 7.90142264483e-9],
    [35, 7, -3.45829826367e-9, 4.71386839716e-9],
    [35, 8, 4.15911228686e-9, 9.21486965423e-9],
    [35, 9, -7.83584593022e-10, -1.08780700595e-9],
    [35, 10, -2.63078124596e-9, 1.14437669825e-8],
    [35, 11, 3.1135284219e-9, -3.11508942142e-9],
    [35, 12, 8.10432165903e-9, -6.4323395678e-9],
    [35, 13, -1.60870380988e-9, 3.02852925442e-9],
    [35, 14, -7.16511186947e-9, -7.02737046917e-9],
    [35, 15, -1.53690564123e-8, 8.75984924717e-9],
    [35, 16, -6.89772047703e-9, -7.36827047584e-9],
    [35, 17, 7.03755899027e-10, -8.82920485773e-9],
    [35, 18, -5.55247661498e-9, -1.14710477959e-8],
    [35, 19, -1.07112499273e-9, -3.41854119412e-9],
    [35, 20, 9.92702305837e-10, -1.13573745208e-10],
    [35, 21, 1.29333785663e-8, -8.17657795386e-10],
    [35, 22, 7.51479477595e-9, 5.7229930908e-9],
    [35, 23, -8.16391242216e-9, -2.22442612532e-9],
    [35, 24, 2.78435090517e-9, 6.38499607176e-9],
    [35, 25, 7.16858934156e-9, 1.99781103645e-9],
    [35, 26, -4.70300232305e-9, 4.61488943108e-9],
    [35, 27, 1.09602089094e-8, -1.33812635796e-8],
    [35, 28, 7.88159460716e-9, -1.53673024839e-8],
    [35, 29, 7.70786810766e-9, 3.40140754669e-9],
    [35, 30, -4.0519283993e-9, 2.87370616224e-9],
    [35, 31, 7.84140204315e-9, 4.0412480788e-9],
    [35, 32, -3.16267901777e-9, -7.41858064221e-9],
    [35, 33, 5.8609633966e-9, -3.07739390905e-9],
    [35, 34, -1.21632099674e-9, 2.66717400938e-9],
    [35, 35, -5.8786572941e-9, -5.01230638002e-9],
    [36, 0, -4.02590604243e-9, 0],
    [36, 1, -1.13386686386e-9, 5.14982653283e-9],
    [36, 2, -4.31575901448e-9, -3.40211031655e-9],
    [36, 3, 7.00409280444e-11, -1.58895672921e-8],
    [36, 4, 3.00961129935e-9, 1.38917218538e-9],
    [36, 5, -7.42261535513e-9, 1.4033786019e-9],
    [36, 6, 1.08546024568e-8, -3.16311943226e-9],
    [36, 7, 1.70813806147e-9, 6.17680210154e-9],
    [36, 8, 3.44939360246e-9, -5.03767857861e-9],
    [36, 9, 2.92192219493e-9, -3.74028113708e-10],
    [36, 10, 4.23119681703e-9, 6.83503143788e-9],
    [36, 11, -4.10039232642e-9, 4.75118294475e-9],
    [36, 12, 4.87204962837e-10, -9.84587714675e-9],
    [36, 13, -6.15416963507e-9, 8.0318113556e-9],
    [36, 14, -1.04141682764e-8, -5.94203574762e-9],
    [36, 15, 9.54892409044e-10, 3.33310574172e-9],
    [36, 16, 1.25505913598e-9, -1.60569406116e-10],
    [36, 17, 4.95066186034e-9, -8.65314022477e-9],
    [36, 18, 1.77184202015e-9, 4.4603340077e-9],
    [36, 19, -5.25149217565e-9, -6.65319486115e-9],
    [36, 20, -6.03793346956e-9, 3.52627660597e-9],
    [36, 21, 1.0690892473e-8, -5.67948915026e-9],
    [36, 22, 3.21356130034e-9, 1.61234121461e-9],
    [36, 23, -3.61160199501e-10, 2.74891917069e-9],
    [36, 24, 2.10662869987e-9, -4.24514998756e-9],
    [36, 25, 4.3497929214e-9, 1.5607147346e-8],
    [36, 26, 3.68762567031e-9, 9.37175113714e-9],
    [36, 27, -7.91229464362e-9, 8.8299681063e-9],
    [36, 28, 2.22637976824e-9, -4.34372617405e-9],
    [36, 29, 1.84511675839e-9, 2.0734471834e-10],
    [36, 30, -1.00411515955e-8, 6.05413293608e-9],
    [36, 31, -8.39084442298e-9, -5.54047445598e-9],
    [36, 32, 1.25654207109e-8, 2.30476235625e-9],
    [36, 33, 3.89957606637e-9, -3.50340856893e-9],
    [36, 34, -9.08693282663e-9, 4.35776976715e-9],
    [36, 35, -1.38812503272e-10, -1.25527291076e-8],
    [36, 36, 4.6014646572e-9, -5.94245336314e-9],
];


}),
"./src/engine/ootk/src/data/values/hpAtmosphere.ts": 
/*!*********************************************************!*\
  !*** ./src/engine/ootk/src/data/values/hpAtmosphere.ts ***!
  \*********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  hpAtmosphere: () => (hpAtmosphere)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Harris-Priester atmosphere data, assuming mean solar flux.
const hpAtmosphere = [
    [100, 4.974e-7, 4.974e-7],
    [120, 2.49e-8, 2.49e-8],
    [130, 8.377e-9, 8.71e-9],
    [140, 3.899e-9, 4.059e-9],
    [150, 2.122e-9, 2.215e-9],
    [160, 1.263e-9, 1.344e-9],
    [170, 8.008e-10, 8.758e-10],
    [180, 5.283e-10, 6.01e-10],
    [190, 3.617e-10, 4.297e-10],
    [200, 2.557e-10, 3.162e-10],
    [210, 1.839e-10, 2.396e-10],
    [220, 1.341e-10, 1.853e-10],
    [230, 9.949e-11, 1.455e-10],
    [240, 7.488e-11, 1.157e-10],
    [250, 5.709e-11, 9.308e-11],
    [260, 4.403e-11, 7.555e-11],
    [270, 3.43e-11, 6.182e-11],
    [280, 2.697e-11, 5.095e-11],
    [290, 2.139e-11, 4.226e-11],
    [300, 1.708e-11, 3.526e-11],
    [320, 1.099e-11, 2.511e-11],
    [340, 7.214e-12, 1.819e-11],
    [360, 4.824e-12, 1.337e-11],
    [380, 3.274e-12, 9.955e-12],
    [400, 2.249e-12, 7.492e-12],
    [420, 1.558e-12, 5.684e-12],
    [440, 1.091e-12, 4.355e-12],
    [460, 7.701e-13, 3.362e-12],
    [480, 5.474e-13, 2.612e-12],
    [500, 3.916e-13, 2.042e-12],
    [520, 2.819e-13, 1.605e-12],
    [540, 2.042e-13, 1.267e-12],
    [560, 1.488e-13, 1.005e-12],
    [580, 1.092e-13, 7.997e-13],
    [600, 8.07e-14, 6.39e-13],
    [620, 6.012e-14, 5.123e-13],
    [640, 4.519e-14, 4.121e-13],
    [660, 3.43e-14, 3.325e-13],
    [680, 2.632e-14, 2.691e-13],
    [700, 2.043e-14, 2.185e-13],
    [720, 1.607e-14, 1.779e-13],
    [740, 1.281e-14, 1.452e-13],
    [760, 1.036e-14, 1.19e-13],
    [780, 8.496e-15, 9.776e-14],
    [800, 7.069e-15, 8.059e-14],
    [840, 4.68e-15, 5.741e-14],
    [880, 3.2e-15, 4.21e-14],
    [920, 2.21e-15, 3.13e-14],
    [960, 1.56e-15, 2.36e-14],
    [1000, 1.15e-15, 1.81e-14],
];


}),
"./src/engine/ootk/src/data/values/iau1980.ts": 
/*!****************************************************!*\
  !*** ./src/engine/ootk/src/data/values/iau1980.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  iau1980: () => (iau1980)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Array of the first 4 IAU-1980 coefficients.
const iau1980 = [
    [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
    [0, 0, 2, -2, 2, -13187, -1.6, 5736, -3.1],
    [0, 0, 2, 0, 2, -2274, -0.2, 977, -0.5],
    [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
];


}),
"./src/engine/ootk/src/enums/AngularDiameterMethod.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/enums/AngularDiameterMethod.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  AngularDiameterMethod: () => (AngularDiameterMethod)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/** Enumeration representing different methods for calculating angular diameter. */
var AngularDiameterMethod;
(function (AngularDiameterMethod) {
    AngularDiameterMethod[AngularDiameterMethod["Circle"] = 0] = "Circle";
    AngularDiameterMethod[AngularDiameterMethod["Sphere"] = 1] = "Sphere";
})(AngularDiameterMethod || (AngularDiameterMethod = {}));


}),
"./src/engine/ootk/src/enums/AngularDistanceMethod.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/enums/AngularDistanceMethod.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  AngularDistanceMethod: () => (AngularDistanceMethod)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/** Enumeration representing different methods for calculating angular distance. */
var AngularDistanceMethod;
(function (AngularDistanceMethod) {
    AngularDistanceMethod[AngularDistanceMethod["Cosine"] = 0] = "Cosine";
    AngularDistanceMethod[AngularDistanceMethod["Haversine"] = 1] = "Haversine";
})(AngularDistanceMethod || (AngularDistanceMethod = {}));


}),
"./src/engine/ootk/src/enums/CatalogSource.ts": 
/*!****************************************************!*\
  !*** ./src/engine/ootk/src/enums/CatalogSource.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CatalogSource: () => (CatalogSource)
});
var CatalogSource;
(function (CatalogSource) {
    CatalogSource["USSF"] = "USSF";
    CatalogSource["CELESTRAK"] = "Celestrak";
    CatalogSource["UNIV_OF_MICH"] = "University of Michigan";
    CatalogSource["CALPOLY"] = "CalPoly";
    CatalogSource["NUSPACE"] = "NuSpace";
    CatalogSource["VIMPEL"] = "JSC Vimpel";
    CatalogSource["TLE_TXT"] = "TLE.txt";
    CatalogSource["EXTRA_JSON"] = "extra.json";
})(CatalogSource || (CatalogSource = {}));


}),
"./src/engine/ootk/src/enums/CommLink.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/enums/CommLink.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CommLink: () => (CommLink)
});
var CommLink;
(function (CommLink) {
    CommLink["AEHF"] = "AEHF";
    CommLink["GALILEO"] = "Galileo";
    CommLink["IRIDIUM"] = "Iridium";
    CommLink["STARLINK"] = "Starlink";
    CommLink["WGS"] = "WGS";
})(CommLink || (CommLink = {}));


}),
"./src/engine/ootk/src/enums/OrbitRegime.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/enums/OrbitRegime.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  OrbitRegime: () => (OrbitRegime)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/** Orbit regime classifications. */
var OrbitRegime;
(function (OrbitRegime) {
    OrbitRegime["LEO"] = "Low Earth Orbit";
    OrbitRegime["MEO"] = "Medium Earth Orbit";
    OrbitRegime["HEO"] = "Highly Eccentric Orbit";
    OrbitRegime["GEO"] = "Geosynchronous Orbit";
    OrbitRegime["OTHER"] = "Uncategorized Orbit";
})(OrbitRegime || (OrbitRegime = {}));


}),
"./src/engine/ootk/src/enums/PassType.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/enums/PassType.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  PassType: () => (PassType)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
var PassType;
(function (PassType) {
    PassType[PassType["OUT_OF_VIEW"] = -1] = "OUT_OF_VIEW";
    PassType[PassType["ENTER"] = 0] = "ENTER";
    PassType[PassType["IN_VIEW"] = 1] = "IN_VIEW";
    PassType[PassType["EXIT"] = 2] = "EXIT";
})(PassType || (PassType = {}));


}),
"./src/engine/ootk/src/enums/Sgp4OpsMode.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/enums/Sgp4OpsMode.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sgp4OpsMode: () => (Sgp4OpsMode)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
var Sgp4OpsMode;
(function (Sgp4OpsMode) {
    Sgp4OpsMode["AFSPC"] = "a";
    Sgp4OpsMode["IMPROVED"] = "i";
})(Sgp4OpsMode || (Sgp4OpsMode = {}));


}),
"./src/engine/ootk/src/enums/index.ts": 
/*!********************************************!*\
  !*** ./src/engine/ootk/src/enums/index.ts ***!
  \********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  AngularDiameterMethod: () => (/* reexport safe */ _AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_2__.AngularDiameterMethod),
  AngularDistanceMethod: () => (/* reexport safe */ _AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_3__.AngularDistanceMethod),
  CatalogSource: () => (/* reexport safe */ _CatalogSource_js__WEBPACK_IMPORTED_MODULE_5__.CatalogSource),
  CommLink: () => (/* reexport safe */ _CommLink_js__WEBPACK_IMPORTED_MODULE_6__.CommLink),
  OrbitRegime: () => (/* reexport safe */ _OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime),
  PassType: () => (/* reexport safe */ _PassType_js__WEBPACK_IMPORTED_MODULE_4__.PassType),
  Sgp4OpsMode: () => (/* reexport safe */ _Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_1__.Sgp4OpsMode)
});
/* ESM import */var _OrbitRegime_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./OrbitRegime.js */ "./src/engine/ootk/src/enums/OrbitRegime.ts");
/* ESM import */var _Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Sgp4OpsMode.js */ "./src/engine/ootk/src/enums/Sgp4OpsMode.ts");
/* ESM import */var _AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./AngularDiameterMethod.js */ "./src/engine/ootk/src/enums/AngularDiameterMethod.ts");
/* ESM import */var _AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./AngularDistanceMethod.js */ "./src/engine/ootk/src/enums/AngularDistanceMethod.ts");
/* ESM import */var _PassType_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./PassType.js */ "./src/engine/ootk/src/enums/PassType.ts");
/* ESM import */var _CatalogSource_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./CatalogSource.js */ "./src/engine/ootk/src/enums/CatalogSource.ts");
/* ESM import */var _CommLink_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./CommLink.js */ "./src/engine/ootk/src/enums/CommLink.ts");









}),
"./src/engine/ootk/src/force/AtmosphericDrag.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/force/AtmosphericDrag.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  AtmosphericDrag: () => (AtmosphericDrag)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Harris-Priester atmospheric drag force model.
 * Atmospheric density model assumes mean solar flux.
 */
class AtmosphericDrag {
    mass;
    area;
    dragCoeff;
    cosine;
    constructor(mass, area, dragCoeff, cosine) {
        this.mass = mass;
        this.area = area;
        this.dragCoeff = dragCoeff;
        this.cosine = cosine;
    }
    static _getHPDensity(state, n) {
        const hpa = _main_js__WEBPACK_IMPORTED_MODULE_0__.DataHandler.getInstance().getHpAtmosphere(state.height);
        if (hpa === null) {
            return 0.0;
        }
        const sunPos = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.positionApparent(state.epoch);
        const sunVec = new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(state.epoch, sunPos, _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin).toITRF().position.normalize();
        const bulVec = sunVec.rotZ(-30.0 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const cosPsi = bulVec.normalize().dot(state.position.normalize());
        const c2Psi2 = 0.5 * (1.0 + cosPsi);
        const cPsi2 = Math.sqrt(c2Psi2);
        const cosPow = cPsi2 > 1e-12 ? c2Psi2 * cPsi2 ** (n - 2) : 0.0;
        const altitude = hpa.height;
        const [h0, min0, max0] = hpa.hp0;
        const [h1, min1, max1] = hpa.hp1;
        const dH = (h0 - altitude) / (h0 - h1);
        const rhoMin = min0 * (min1 / min0) ** dH;
        if (cosPow === 0) {
            return rhoMin;
        }
        const rhoMax = max0 * (max1 / max0) ** dH;
        return rhoMin + (rhoMax - rhoMin) * cosPow;
    }
    acceleration(state) {
        const itrfState = state.toITRF();
        const density = AtmosphericDrag._getHPDensity(itrfState, this.cosine);
        if (density === 0) {
            return _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin;
        }
        const rotation = new _main_js__WEBPACK_IMPORTED_MODULE_0__.ITRF(state.epoch, _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.rotation, _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin).toJ2000().position;
        const vRel = state.velocity.subtract(rotation.cross(state.position))
            .scale(1000.0);
        const vm = vRel.magnitude();
        const fScale = -0.5 * density * ((this.dragCoeff * this.area) / this.mass) * vm;
        return vRel.scale(fScale / 1000.0);
    }
}


}),
"./src/engine/ootk/src/force/EarthGravity.ts": 
/*!***************************************************!*\
  !*** ./src/engine/ootk/src/force/EarthGravity.ts ***!
  \***************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EarthGravity: () => (EarthGravity)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * designed to model the Earth's gravitational field, which is not uniformly distributed due to variations in mass
 * distribution within the Earth and the Earth's shape (it's not a perfect sphere). To accurately model this complex
 * field, the gravity model is expanded into a series of spherical harmonics, characterized by their degree and order.
 *
 * This `degree` parameter is related to the spatial resolution of the gravity model. A higher degree corresponds to a
 * finer resolution, capable of representing smaller-scale variations in the gravity field. The degree essentially
 * denotes how many times the gravitational potential function varies over the surface of the Earth.
 *
 * For each degree, there can be multiple orders ranging from 0 up to the degree. The `order` accounts for the
 * longitudinal variation in the gravity field. Each order within a degree captures different characteristics of the
 * gravity anomalies.
 *
 * `Degree 0` corresponds to the overall, mean gravitational force of the Earth (considered as a point mass).
 *
 * `Degree 1` terms are related to the Earth's center of mass but are usually not used because the center of mass is
 * defined as the origin of the coordinate system.
 *
 * `Degree 2` and higher capture the deviations from this spherical symmetry, such as the flattening at the poles and
 * bulging at the equator (degree 2), and other anomalies at finer scales as the degree increases.
 */
class EarthGravity {
    degree;
    order;
    _asphericalFlag;
    /**
     * Creates a new instance of the EarthGravity class.
     * @param degree The degree of the Earth's gravity field. Must be between 0 and 36.
     * @param order The order of the Earth's gravity field. Must be between 0 and 36.
     */
    constructor(degree, order) {
        this.degree = Math.min(Math.max(degree, 0), 36);
        this.order = Math.min(Math.max(order, 0), 36);
        this._asphericalFlag = degree >= 2;
    }
    _spherical(state) {
        const rMag = state.position.magnitude();
        return state.position.scale(-_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / (rMag * rMag * rMag));
    }
    // eslint-disable-next-line max-statements
    _aspherical(state) {
        const posEcef = state.toITRF().position;
        const ri = 1.0 / posEcef.magnitude();
        const xor = posEcef.x * ri;
        const yor = posEcef.y * ri;
        const zor = posEcef.z * ri;
        const ep = zor;
        const reor = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusEquator * ri;
        let reorn = reor;
        const muor2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu * ri * ri;
        let sumH = 0.0;
        let sumGm = 0.0;
        let sumJ = 0.0;
        let sumK = 0.0;
        const cTil = new Float64Array(this.order + 4);
        const sTil = new Float64Array(this.order + 4);
        const pN = new Float64Array(this.order + 4);
        const pNm1 = new Float64Array(this.order + 4);
        const pNm2 = new Float64Array(this.order + 4);
        pNm2[0] = 1.0;
        pNm1[0] = ep;
        pNm1[1] = 1.0;
        cTil[0] = 1.0;
        cTil[1] = xor;
        sTil[1] = yor;
        const dh = _main_js__WEBPACK_IMPORTED_MODULE_0__.DataHandler.getInstance();
        for (let n = 2, nm1 = 1, nm2 = 0, np1 = 3; n <= this.degree; nm2++, nm1++, n++, np1++) {
            const twonm1 = 2.0 * n - 1.0;
            reorn *= reor;
            const cN0 = dh.getEgm96Coeffs(n, 0)[2];
            pN[0] = (twonm1 * ep * pNm1[0] - nm1 * pNm2[0]) / n;
            pN[1] = pNm2[1] + twonm1 * pNm1[0];
            pN[2] = pNm2[2] + twonm1 * pNm1[1];
            let sumHn = pN[1] * cN0;
            let sumGmn = pN[0] * cN0 * np1;
            if (this.order > 0) {
                let sumJn = 0.0;
                let sumKn = 0.0;
                cTil[n] = cTil[1] * cTil[nm1] - sTil[1] * sTil[nm1];
                sTil[n] = sTil[1] * cTil[nm1] + cTil[1] * sTil[nm1];
                const lim = n < this.order ? n : this.order;
                for (let m = 1, mm1 = 0, mm2 = -1, mp1 = 2, mp2 = 3; m <= lim; mm2++, mm1++, m++, mp1++, mp2++) {
                    pN[mp1] = pNm2[mp1] + twonm1 * pNm1[m];
                    const dm = m;
                    const npmp1 = n + mp1;
                    const pNm = pN[m];
                    const pNmp1 = pN[mp1];
                    const coefs = dh.getEgm96Coeffs(n, m);
                    const cNm = coefs[2];
                    const sNm = coefs[3];
                    const mxPnm = dm * pNm;
                    const bNmtil = cNm * cTil[m] + sNm * sTil[m];
                    const pNmBnm = pNm * bNmtil;
                    const bNmtm1 = cNm * cTil[mm1] + sNm * sTil[mm1];
                    const aNmtm1 = cNm * sTil[mm1] - sNm * cTil[mm1];
                    sumHn += pNmp1 * bNmtil;
                    sumGmn += npmp1 * pNmBnm;
                    sumJn += mxPnm * bNmtm1;
                    sumKn -= mxPnm * aNmtm1;
                }
                sumJ += reorn * sumJn;
                sumK += reorn * sumKn;
            }
            sumH += reorn * sumHn;
            sumGm += reorn * sumGmn;
            if (n < this.degree) {
                for (let i = 0; i <= n; i++) {
                    pNm2[i] = pNm1[i];
                    pNm1[i] = pN[i];
                }
            }
        }
        const lambda = sumGm + ep * sumH;
        const g = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(-muor2 * (lambda * xor - sumJ), -muor2 * (lambda * yor - sumK), -muor2 * (lambda * zor - sumH));
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.ITRF(state.epoch, g, _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin).toJ2000().position;
    }
    acceleration(state) {
        let accVec = this._spherical(state);
        if (this._asphericalFlag) {
            accVec = accVec.add(this._aspherical(state));
        }
        return accVec;
    }
}


}),
"./src/engine/ootk/src/force/Force.ts": 
/*!********************************************!*\
  !*** ./src/engine/ootk/src/force/Force.ts ***!
  \********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Force: () => (Force)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable @typescript-eslint/no-unused-vars */
/* eslint-disable class-methods-use-this */
// / Base class for perturbation forces.
class Force {
    /**
     * Calculate the acceleration due to the perturbing force on a given
     * state vector.
     * @param state The state vector.
     * @throws If the force cannot be calculated.
     */
    acceleration(state) {
        throw Error('Not implemented');
    }
}


}),
"./src/engine/ootk/src/force/ForceModel.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/force/ForceModel.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ForceModel: () => (ForceModel)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _AtmosphericDrag_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./AtmosphericDrag.js */ "./src/engine/ootk/src/force/AtmosphericDrag.ts");
/* ESM import */var _EarthGravity_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./EarthGravity.js */ "./src/engine/ootk/src/force/EarthGravity.ts");
/* ESM import */var _Gravity_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Gravity.js */ "./src/engine/ootk/src/force/Gravity.ts");
/* ESM import */var _SolarRadiationPressure_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./SolarRadiationPressure.js */ "./src/engine/ootk/src/force/SolarRadiationPressure.ts");
/* ESM import */var _ThirdBodyGravity_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./ThirdBodyGravity.js */ "./src/engine/ootk/src/force/ThirdBodyGravity.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






// / Force model for spacecraft propagation.
class ForceModel {
    _centralGravity;
    _thirdBodyGravity;
    _solarRadiationPressure;
    _atmosphericDrag;
    _maneuverThrust = null;
    setGravity(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        this._centralGravity = new _Gravity_js__WEBPACK_IMPORTED_MODULE_3__.Gravity(mu);
        return this;
    }
    setEarthGravity(degree, order) {
        this._centralGravity = new _EarthGravity_js__WEBPACK_IMPORTED_MODULE_2__.EarthGravity(degree, order);
    }
    setThirdBodyGravity({ moon = false, sun = false }) {
        this._thirdBodyGravity = new _ThirdBodyGravity_js__WEBPACK_IMPORTED_MODULE_5__.ThirdBodyGravity(moon, sun);
    }
    setSolarRadiationPressure(mass, area, coeff = 1.2) {
        this._solarRadiationPressure = new _SolarRadiationPressure_js__WEBPACK_IMPORTED_MODULE_4__.SolarRadiationPressure(mass, area, coeff);
    }
    /**
     * Sets the atmospheric drag for the force model.
     * @deprecated This is still a work in progress!
     * @param mass - The mass of the object.
     * @param area - The cross-sectional area of the object.
     * @param coeff - The drag coefficient. Default value is 2.2.
     * @param cosine - The cosine of the angle between the object's velocity vector and the drag force vector.
     */
    setAtmosphericDrag(mass, area, coeff = 2.2, cosine = 4) {
        this._atmosphericDrag = new _AtmosphericDrag_js__WEBPACK_IMPORTED_MODULE_1__.AtmosphericDrag(mass, area, coeff, cosine);
    }
    loadManeuver(maneuver) {
        this._maneuverThrust = maneuver;
    }
    clearManeuver() {
        this._maneuverThrust = null;
    }
    acceleration(state) {
        let accVec = _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin;
        if (this._centralGravity) {
            accVec = accVec.add(this._centralGravity.acceleration(state));
        }
        if (this._thirdBodyGravity) {
            accVec = accVec.add(this._thirdBodyGravity.acceleration(state));
        }
        if (this._solarRadiationPressure) {
            accVec = accVec.add(this._solarRadiationPressure.acceleration(state));
        }
        if (this._atmosphericDrag) {
            accVec = accVec.add(this._atmosphericDrag.acceleration(state));
        }
        if (this._maneuverThrust) {
            accVec = accVec.add(this._maneuverThrust.acceleration(state));
        }
        return accVec;
    }
    derivative(state) {
        return state.velocity.join(this.acceleration(state));
    }
}


}),
"./src/engine/ootk/src/force/Gravity.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/force/Gravity.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Gravity: () => (Gravity)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Simple central-body gravity model.
class Gravity {
    // / Gravitational parameter _(km²/s³)_.
    mu;
    /**
     * Create a new [Gravity] object with optional gravitational
     * @param mu Gravitational parameter _(km²/s³)_.
     */
    constructor(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        this.mu = mu;
    }
    /**
     * Calculates the gravitational force in spherical coordinates.
     * @param state The J2000 state containing the position and velocity vectors.
     * @returns The gravitational force vector in spherical coordinates.
     */
    _spherical(state) {
        const rMag = state.position.magnitude();
        return state.position.scale(-this.mu / (rMag * rMag * rMag));
    }
    /**
     * Calculates the acceleration due to gravity at a given state.
     * @param state The J2000 state at which to calculate the acceleration.
     * @returns The acceleration vector.
     */
    acceleration(state) {
        return this._spherical(state);
    }
}


}),
"./src/engine/ootk/src/force/SolarRadiationPressure.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/force/SolarRadiationPressure.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  SolarRadiationPressure: () => (SolarRadiationPressure)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _Force_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Force.js */ "./src/engine/ootk/src/force/Force.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


// / Solar radiation pressure model.
class SolarRadiationPressure extends _Force_js__WEBPACK_IMPORTED_MODULE_1__.Force {
    mass;
    area;
    reflectCoeff;
    // / Create a new [SolarRadiationPressure] object.
    constructor(mass, area, reflectCoeff) {
        super();
        this.mass = mass;
        this.area = area;
        this.reflectCoeff = reflectCoeff;
    }
    // / Solar pressure _(N/m²)_;
    static _kRef = 4.56e-6 * _main_js__WEBPACK_IMPORTED_MODULE_0__.astronomicalUnit ** 2;
    acceleration(state) {
        const rSun = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.positionApparent(state.epoch);
        const r = state.position.subtract(rSun);
        const rMag = r.magnitude();
        const r2 = rMag * rMag;
        const ratio = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.lightingRatio(state.position, rSun);
        const p = (ratio * SolarRadiationPressure._kRef) / r2;
        const flux = r.scale(p / rMag);
        return flux.scale(((this.area * this.reflectCoeff) / this.mass) * 1e-3);
    }
}


}),
"./src/engine/ootk/src/force/ThirdBodyGravity.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/force/ThirdBodyGravity.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ThirdBodyGravity: () => (ThirdBodyGravity)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Third-body gravity model.
class ThirdBodyGravity {
    moon;
    sun;
    // / Create a new [ThirdBodyGravity] object with the selected bodies enabled.
    constructor(moon = false, sun = false) {
        this.moon = moon;
        this.sun = sun;
        // Nothing to do here.
    }
    static _moonGravity(state) {
        const rMoon = _main_js__WEBPACK_IMPORTED_MODULE_0__.Moon.eci(state.epoch);
        const aNum = rMoon.subtract(state.position);
        const aDen = aNum.magnitude() ** 3;
        const bNum = rMoon;
        const bDen = rMoon.magnitude() ** 3;
        const gravity = aNum.scale(1 / aDen).add(bNum.scale(-1 / bDen));
        return gravity.scale(_main_js__WEBPACK_IMPORTED_MODULE_0__.Moon.mu);
    }
    static _sunGravity(state) {
        const rSun = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.positionApparent(state.epoch);
        const aNum = rSun.subtract(state.position);
        const aDen = aNum.magnitude() ** 3;
        const bNum = rSun;
        const bDen = rSun.magnitude() ** 3;
        const gravity = aNum.scale(1 / aDen).add(bNum.scale(-1 / bDen));
        return gravity.scale(_main_js__WEBPACK_IMPORTED_MODULE_0__.Sun.mu);
    }
    acceleration(state) {
        let accVec = _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin;
        if (this.moon) {
            accVec = accVec.add(ThirdBodyGravity._moonGravity(state));
        }
        if (this.sun) {
            accVec = accVec.add(ThirdBodyGravity._sunGravity(state));
        }
        return accVec;
    }
}


}),
"./src/engine/ootk/src/force/Thrust.ts": 
/*!*********************************************!*\
  !*** ./src/engine/ootk/src/force/Thrust.ts ***!
  \*********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Thrust: () => (Thrust)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Thrust force model.
class Thrust {
    center;
    radial;
    intrack;
    crosstrack;
    durationRate;
    constructor(center, radial, intrack, crosstrack, durationRate = 0.0) {
        this.center = center;
        this.radial = radial;
        this.intrack = intrack;
        this.crosstrack = crosstrack;
        this.durationRate = durationRate;
        this.deltaV = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(radial * 1e-3, intrack * 1e-3, crosstrack * 1e-3);
    }
    deltaV;
    get magnitude() {
        return this.deltaV.magnitude() * 1000.0;
    }
    get duration() {
        return this.magnitude * this.durationRate;
    }
    get start() {
        return this.center.roll(-0.5 * this.duration);
    }
    get stop() {
        return this.center.roll(0.5 * this.duration);
    }
    acceleration(state) {
        const relative = new _main_js__WEBPACK_IMPORTED_MODULE_0__.RIC(_main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin, this.deltaV.scale(1.0 / this.duration));
        return relative.toJ2000(state).velocity.subtract(state.velocity);
    }
    apply(state) {
        const relative = new _main_js__WEBPACK_IMPORTED_MODULE_0__.RIC(_main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D.origin, this.deltaV);
        return relative.toJ2000(state);
    }
    get isImpulsive() {
        return this.duration <= 0;
    }
}


}),
"./src/engine/ootk/src/force/index.ts": 
/*!********************************************!*\
  !*** ./src/engine/ootk/src/force/index.ts ***!
  \********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ForceModel: () => (/* reexport safe */ _ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel)
});
/* ESM import */var _ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");



}),
"./src/engine/ootk/src/interfaces/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/interfaces/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);



}),
"./src/engine/ootk/src/interpolator/CubicSpline.ts": 
/*!*********************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/CubicSpline.ts ***!
  \*********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CubicSpline: () => (CubicSpline)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Container for cubic spline data.
class CubicSpline {
    t0;
    p0;
    m0;
    t1;
    p1;
    m1;
    // / Create a new [CubicSpline] object.
    constructor(t0, p0, m0, t1, p1, m1) {
        this.t0 = t0;
        this.p0 = p0;
        this.m0 = m0;
        this.t1 = t1;
        this.p1 = p1;
        this.m1 = m1;
        // Nothing to do here.
    }
    // / Interpolate position at the provided time [t] _(POSIX seconds)_.
    position_(t) {
        const t2 = t * t;
        const t3 = t2 * t;
        const r0 = this.p0.scale(2 * t3 - 3 * t2 + 1);
        const v0 = this.m0.scale((t3 - 2 * t2 + t) * (this.t1 - this.t0));
        const r1 = this.p1.scale(-2 * t3 + 3 * t2);
        const v1 = this.m1.scale((t3 - t2) * (this.t1 - this.t0));
        return r0.add(v0).add(r1).add(v1);
    }
    // / Interpolate velocity at the provided time [t] _(POSIX seconds)_.
    velocity_(t) {
        const t2 = t * t;
        const r0 = this.p0.scale(6 * t2 - 6 * t);
        const v0 = this.m0.scale((3 * t2 - 4 * t + 1) * (this.t1 - this.t0));
        const r1 = this.p1.scale(-6 * t2 + 6 * t);
        const v1 = this.m1.scale((3 * t2 - 2 * t) * (this.t1 - this.t0));
        return r0
            .add(v0)
            .add(r1)
            .add(v1)
            .scale(1 / (this.t1 - this.t0));
    }
    /**
     * Interpolates the position and velocity at a given time.
     * (km) and velocity (km/s) vectors at the provided time.
     * @param t The time value to interpolate at _(POSIX seconds)_.
     * @returns An array containing the interpolated position and velocity as Vector3D objects.
     */
    interpolate(t) {
        const n = (t - this.t0) / (this.t1 - this.t0);
        return [this.position_(n), this.velocity_(n)];
    }
}


}),
"./src/engine/ootk/src/interpolator/CubicSplineInterpolator.ts": 
/*!*********************************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/CubicSplineInterpolator.ts ***!
  \*********************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  CubicSplineInterpolator: () => (CubicSplineInterpolator)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _CubicSpline_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./CubicSpline.js */ "./src/engine/ootk/src/interpolator/CubicSpline.ts");
/* ESM import */var _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./StateInterpolator.js */ "./src/engine/ootk/src/interpolator/StateInterpolator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



/**
 * Cubic spline ephemeris interpolator.
 *
 * The [CubicSplineInterpolator] is a very fast and accurate interpolator
 * at the expense of memory due to the cached spline pairs used in the
 * interpolation operation. Accuracy is significantly impacted when using
 * sparse ephemerides.
 */
class CubicSplineInterpolator extends _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_2__.StateInterpolator {
    splines_;
    constructor(splines_) {
        super();
        this.splines_ = splines_;
    }
    static fromEphemeris(ephemeris) {
        const splines = [];
        for (let i = 0; i < ephemeris.length - 1; i++) {
            const e0 = ephemeris[i];
            const t0 = e0.epoch.posix;
            const p0 = e0.position;
            const m0 = e0.velocity;
            const e1 = ephemeris[i + 1];
            const t1 = e1.epoch.posix;
            const p1 = e1.position;
            const m1 = e1.velocity;
            splines.push(new _CubicSpline_js__WEBPACK_IMPORTED_MODULE_1__.CubicSpline(t0, p0, m0, t1, p1, m1));
        }
        return new CubicSplineInterpolator(splines);
    }
    get sizeBytes() {
        return (64 * 14 * this.splines_.length) / 8;
    }
    matchSpline_(posix) {
        let left = 0;
        let right = this.splines_.length;
        while (left < right) {
            const middle = (left + right) >> 1;
            if (this.splines_[middle].t1 < posix) {
                left = middle + 1;
            }
            else {
                right = middle;
            }
        }
        return this.splines_[left];
    }
    interpolate(epoch) {
        if (!this.inWindow(epoch)) {
            return null;
        }
        const posix = epoch.posix;
        const splineVecs = this
            .matchSpline_(posix)
            .interpolate(posix);
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(epoch, splineVecs[0], splineVecs[1]);
    }
    window() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochWindow(new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(this.splines_[0].t0), new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(this.splines_[this.splines_.length - 1].t1));
    }
}


}),
"./src/engine/ootk/src/interpolator/Interpolator.ts": 
/*!**********************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/Interpolator.ts ***!
  \**********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Interpolator: () => (Interpolator)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Interpolator base class.
class Interpolator {
    /*
     * Return `true` if the provided [epoch] is within this interpolator's
     * cached value range.
     */
    inWindow(epoch) {
        const start = this.window().start;
        const stop = this.window().end;
        return start <= epoch && epoch <= stop;
    }
    /*
     * Calculate the start/stop epoch between this and another [Interpolator].
     *
     * Returns `null` if there is no overlap between interpolators.
     */
    overlap(interpolator) {
        const x1 = this.window().start;
        const x2 = this.window().end;
        const y1 = interpolator.window().start;
        const y2 = interpolator.window().end;
        if (x1 <= y2 && y1 <= x2) {
            const e1 = new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(Math.max(x1.posix, y1.posix));
            const e2 = new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(Math.min(x2.posix, y2.posix));
            return new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochWindow(e1, e2);
        }
        return null;
    }
}


}),
"./src/engine/ootk/src/interpolator/LagrangeInterpolator.ts": 
/*!******************************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/LagrangeInterpolator.ts ***!
  \******************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  LagrangeInterpolator: () => (LagrangeInterpolator)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./StateInterpolator.js */ "./src/engine/ootk/src/interpolator/StateInterpolator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


class LagrangeInterpolator extends _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_1__.StateInterpolator {
    t_;
    x_;
    y_;
    z_;
    order;
    constructor(t, x, y, z, order = 10) {
        super();
        this.t_ = t;
        this.x_ = x;
        this.y_ = y;
        this.z_ = z;
        this.order = order;
    }
    /**
     * Creates a LagrangeInterpolator from an array of J2000 ephemeris data.
     * @param ephemeris - The array of J2000 ephemeris data.
     * @param order - The order of the LagrangeInterpolator. Default is 10.
     * @returns A new LagrangeInterpolator instance.
     */
    static fromEphemeris(ephemeris, order = 10) {
        const k = ephemeris.length;
        const t = new Float64Array(k);
        const x = new Float64Array(k);
        const y = new Float64Array(k);
        const z = new Float64Array(k);
        for (let i = 0; i < k; i++) {
            const state = ephemeris[i];
            t[i] = state.epoch.posix;
            x[i] = state.position.x;
            y[i] = state.position.y;
            z[i] = state.position.z;
        }
        return new LagrangeInterpolator(t, x, y, z, order);
    }
    get sizeBytes() {
        return (64 * 4 * this.t_.length) / 8;
    }
    interpolate(epoch) {
        if (!this.inWindow(epoch)) {
            return null;
        }
        const posix = epoch.posix;
        const subDex = this.slice_(posix);
        const start = subDex.left;
        const stop = subDex.right;
        const ts = this.t_.subarray(start, stop);
        const xs = this.x_.subarray(start, stop);
        const ys = this.y_.subarray(start, stop);
        const zs = this.z_.subarray(start, stop);
        const position = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(LagrangeInterpolator.position_(ts, xs, posix), LagrangeInterpolator.position_(ts, ys, posix), LagrangeInterpolator.position_(ts, zs, posix));
        const velocity = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(LagrangeInterpolator.velocity_(ts, xs, posix), LagrangeInterpolator.velocity_(ts, ys, posix), LagrangeInterpolator.velocity_(ts, zs, posix));
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(epoch, position, velocity);
    }
    static position_(xs, ys, x) {
        const k = xs.length - 1;
        let result = 0.0;
        for (let j = 0; j < k; j++) {
            let product = ys[j];
            for (let m = 0; m < k; m++) {
                if (j === m) {
                    continue;
                }
                product *= (x - xs[m]) / (xs[j] - xs[m]);
            }
            result += product;
        }
        return result;
    }
    static velocity_(xs, ys, x) {
        const k = xs.length;
        let result = 0.0;
        for (let j = 0; j < k; j++) {
            let total = 0.0;
            for (let i = 0; i < k; i++) {
                if (i === j) {
                    continue;
                }
                let product = 1 / (xs[j] - xs[i]);
                for (let m = 0; m < k; m++) {
                    if (m === i || m === j) {
                        continue;
                    }
                    product *= (x - xs[m]) / (xs[j] - xs[m]);
                }
                total += product;
            }
            result += ys[j] * total;
        }
        return result;
    }
    static _getClosest(target, t1, d1, t2, d2) {
        return target - t1 >= t2 - target ? d2 : d1;
    }
    slice_(posix) {
        const n = this.t_.length;
        if (posix <= this.t_[0]) {
            return { left: 0, right: this.order };
        }
        if (posix >= this.t_[n - 1]) {
            return { left: n - this.order, right: n };
        }
        let i = 0;
        let j = this.t_.length;
        let mid = 0;
        while (i < j) {
            mid = (i + j) >> 1;
            if (this.t_[mid] === posix) {
                break;
            }
            if (posix < this.t_[mid]) {
                if (mid > 0 && posix > this.t_[mid - 1]) {
                    mid = LagrangeInterpolator._getClosest(posix, this.t_[mid - 1], mid - 1, this.t_[mid], mid);
                    break;
                }
                j = mid;
            }
            else {
                if (mid < this.t_.length - 1 && posix < this.t_[mid + 1]) {
                    mid = LagrangeInterpolator._getClosest(posix, this.t_[mid], mid, this.t_[mid + 1], mid + 1);
                    break;
                }
                i = mid + 1;
            }
        }
        const offset = Math.floor(this.order / 2);
        const left = mid - offset;
        const right = mid + offset - (this.order % 2 === 1 ? 1 : 0);
        if (left < 0) {
            return { left: 0, right: this.order };
        }
        if (right > n) {
            return { left: n - this.order, right: n };
        }
        return { left, right };
    }
    window() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochWindow(new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(this.t_[0]), new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochUTC(this.t_[this.t_.length - 1]));
    }
}


}),
"./src/engine/ootk/src/interpolator/StateInterpolator.ts": 
/*!***************************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/StateInterpolator.ts ***!
  \***************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  StateInterpolator: () => (StateInterpolator)
});
/* ESM import */var _Interpolator_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Interpolator.js */ "./src/engine/ootk/src/interpolator/Interpolator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable @typescript-eslint/no-unused-vars */
/* eslint-disable class-methods-use-this */

// / Base class for state vector interpolators.
class StateInterpolator extends _Interpolator_js__WEBPACK_IMPORTED_MODULE_0__.Interpolator {
    /**
     * Interpolates the state at the given epoch.
     * @param epoch The epoch in UTC format.
     * @throws If the interpolator has not been initialized.
     */
    interpolate(epoch) {
        throw new Error('Not implemented.');
    }
    // / Return the size _(bytes)_ of this interpolator's cached data.
    get sizeBytes() {
        throw new Error('Not implemented.');
    }
}


}),
"./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts": 
/*!*********************************************************************!*\
  !*** ./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts ***!
  \*********************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  VerletBlendInterpolator: () => (VerletBlendInterpolator)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _CubicSplineInterpolator_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./CubicSplineInterpolator.js */ "./src/engine/ootk/src/interpolator/CubicSplineInterpolator.ts");
/* ESM import */var _LagrangeInterpolator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./LagrangeInterpolator.js */ "./src/engine/ootk/src/interpolator/LagrangeInterpolator.ts");
/* ESM import */var _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./StateInterpolator.js */ "./src/engine/ootk/src/interpolator/StateInterpolator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




/**
 * Two-body Velocity Verlet Blend interpolator.
 *
 * The [VerletBlendInterpolator] retains the original ephemerides, so the
 * original _"truth"_ states can be retrieved if needed without imparting any
 * additional error, so this can be used to build other interpolator types.
 * The implementation is simple and very tolerant when working with sparse
 * ephemerides.
 */
class VerletBlendInterpolator extends _StateInterpolator_js__WEBPACK_IMPORTED_MODULE_3__.StateInterpolator {
    ephemeris;
    constructor(ephemeris) {
        super();
        this.ephemeris = ephemeris;
    }
    get sizeBytes() {
        return (64 * 7 * this.ephemeris.length) / 8;
    }
    window() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.EpochWindow(this.ephemeris[0].epoch, this.ephemeris[this.ephemeris.length - 1].epoch);
    }
    static getClosest_(target, s1, s2) {
        return target - s1.epoch.posix >= s2.epoch.posix - target ? s2 : s1;
    }
    matchState_(epoch) {
        const target = epoch.posix;
        if (target <= this.ephemeris[0].epoch.posix) {
            return this.ephemeris[0];
        }
        if (target >= this.ephemeris[this.ephemeris.length - 1].epoch.posix) {
            return this.ephemeris[this.ephemeris.length - 1];
        }
        let i = 0;
        let j = this.ephemeris.length;
        let mid = 0;
        while (i < j) {
            mid = (i + j) >> 1;
            if (this.ephemeris[mid].epoch.posix === target) {
                return this.ephemeris[mid];
            }
            if (target < this.ephemeris[mid].epoch.posix) {
                if (mid > 0 && target > this.ephemeris[mid - 1].epoch.posix) {
                    return VerletBlendInterpolator.getClosest_(target, this.ephemeris[mid - 1], this.ephemeris[mid]);
                }
                j = mid;
            }
            else {
                if (mid < this.ephemeris.length - 1 && target < this.ephemeris[mid + 1].epoch.posix) {
                    return VerletBlendInterpolator.getClosest_(target, this.ephemeris[mid], this.ephemeris[mid + 1]);
                }
                i = mid + 1;
            }
        }
        return this.ephemeris[mid];
    }
    static _gravity(position) {
        const r = position.magnitude();
        return position.scale(-_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu / (r * r * r));
    }
    static integrate_(state, step) {
        const x0 = state.position;
        const a0 = VerletBlendInterpolator._gravity(x0);
        const v0 = state.velocity;
        const x1 = x0
            .add(v0.scale(step))
            .add(a0.scale(0.5 * step * step));
        const a1 = VerletBlendInterpolator._gravity(x1);
        const v1 = v0
            .add(a0.add(a1).scale(0.5 * step));
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(state.epoch.roll(step), x1, v1);
    }
    interpolate(epoch) {
        if (!this.inWindow(epoch)) {
            return null;
        }
        let state = this.matchState_(epoch);
        while (state.epoch.posix !== epoch.posix) {
            const delta = epoch.posix - state.epoch.posix;
            const stepMag = Math.min(5.0, Math.abs(delta));
            const stepSize = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.copySign)(stepMag, delta);
            state = VerletBlendInterpolator.integrate_(state, stepSize);
        }
        return state;
    }
    getCachedState(epoch) {
        if (!this.inWindow(epoch)) {
            return null;
        }
        return this.matchState_(epoch);
    }
    toCubicSpline() {
        return _CubicSplineInterpolator_js__WEBPACK_IMPORTED_MODULE_1__.CubicSplineInterpolator.fromEphemeris(this.ephemeris);
    }
    toLagrange(order = 10) {
        return _LagrangeInterpolator_js__WEBPACK_IMPORTED_MODULE_2__.LagrangeInterpolator.fromEphemeris(this.ephemeris, order);
    }
}


}),
"./src/engine/ootk/src/main.ts": 
/*!*************************************!*\
  !*** ./src/engine/ootk/src/main.ts ***!
  \*************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  AngularDiameterMethod: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod),
  AngularDistanceMethod: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.AngularDistanceMethod),
  BaseObject: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.BaseObject),
  BatchLeastSquaresOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.BatchLeastSquaresOD),
  BatchLeastSquaresResult: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.BatchLeastSquaresResult),
  BoxMuller: () => (/* reexport safe */ _operations_index_js__WEBPACK_IMPORTED_MODULE_13__.BoxMuller),
  CatalogSource: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.CatalogSource),
  Celestial: () => (/* reexport safe */ _body_index_js__WEBPACK_IMPORTED_MODULE_8__.Celestial),
  ClassicalElements: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.ClassicalElements),
  CommLink: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.CommLink),
  CovarianceFrame: () => (/* reexport safe */ _covariance_index_js__WEBPACK_IMPORTED_MODULE_17__.CovarianceFrame),
  CovarianceSample: () => (/* reexport safe */ _covariance_index_js__WEBPACK_IMPORTED_MODULE_17__.CovarianceSample),
  DEG2RAD: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.DEG2RAD),
  DataHandler: () => (/* reexport safe */ _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_11__.DataHandler),
  DetailedSatellite: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.DetailedSatellite),
  DetailedSensor: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.DetailedSensor),
  Earth: () => (/* reexport safe */ _body_index_js__WEBPACK_IMPORTED_MODULE_8__.Earth),
  Epoch: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.Epoch),
  EpochGPS: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochGPS),
  EpochTAI: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochTAI),
  EpochTDB: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochTDB),
  EpochTT: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochTT),
  EpochUTC: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochUTC),
  EpochWindow: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.EpochWindow),
  EquinoctialElements: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.EquinoctialElements),
  EulerAngles: () => (/* reexport safe */ _operations_index_js__WEBPACK_IMPORTED_MODULE_13__.EulerAngles),
  ForceModel: () => (/* reexport safe */ _force_index_js__WEBPACK_IMPORTED_MODULE_14__.ForceModel),
  FormatTle: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.FormatTle),
  Geodetic: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.Geodetic),
  GibbsIOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.GibbsIOD),
  GoodingIOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.GoodingIOD),
  GroundObject: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.GroundObject),
  HerrickGibbsIOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.HerrickGibbsIOD),
  Hill: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.Hill),
  ITRF: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.ITRF),
  J2000: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.J2000),
  LambertIOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.LambertIOD),
  LandObject: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.LandObject),
  MILLISECONDS_PER_DAY: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.MILLISECONDS_PER_DAY),
  MILLISECONDS_PER_SECOND: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.MILLISECONDS_PER_SECOND),
  MILLISECONDS_TO_DAYS: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.MILLISECONDS_TO_DAYS),
  MINUTES_PER_DAY: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.MINUTES_PER_DAY),
  MS_PER_DAY: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.MS_PER_DAY),
  Marker: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.Marker),
  Matrix: () => (/* reexport safe */ _operations_operations_js__WEBPACK_IMPORTED_MODULE_6__.Matrix),
  ModifiedGoodingIOD: () => (/* reexport safe */ _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__.ModifiedGoodingIOD),
  Moon: () => (/* reexport safe */ _body_index_js__WEBPACK_IMPORTED_MODULE_8__.Moon),
  OrbitRegime: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.OrbitRegime),
  PI: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.PI),
  PassType: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.PassType),
  PayloadStatus: () => (/* reexport safe */ _types_types_js__WEBPACK_IMPORTED_MODULE_1__.PayloadStatus),
  Quaternion: () => (/* reexport safe */ _operations_index_js__WEBPACK_IMPORTED_MODULE_13__.Quaternion),
  RAD2DEG: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.RAD2DEG),
  RADIUS_OF_EARTH: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.RADIUS_OF_EARTH),
  RAE: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.RAE),
  RIC: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.RIC),
  RadecGeocentric: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.RadecGeocentric),
  RadecTopocentric: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.RadecTopocentric),
  Random: () => (/* reexport safe */ _operations_operations_js__WEBPACK_IMPORTED_MODULE_6__.Random),
  RandomGaussianSource: () => (/* reexport safe */ _operations_index_js__WEBPACK_IMPORTED_MODULE_13__.RandomGaussianSource),
  RelativeState: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.RelativeState),
  RfSensor: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.RfSensor),
  RungeKutta4Propagator: () => (/* reexport safe */ _propagator_index_js__WEBPACK_IMPORTED_MODULE_15__.RungeKutta4Propagator),
  RungeKutta89Propagator: () => (/* reexport safe */ _propagator_index_js__WEBPACK_IMPORTED_MODULE_15__.RungeKutta89Propagator),
  Satellite: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.Satellite),
  Sensor: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.Sensor),
  Sgp4: () => (/* reexport safe */ _sgp4_index_js__WEBPACK_IMPORTED_MODULE_12__.Sgp4),
  Sgp4OpsMode: () => (/* reexport safe */ _enums_index_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode),
  Sgp4Propagator: () => (/* reexport safe */ _propagator_index_js__WEBPACK_IMPORTED_MODULE_15__.Sgp4Propagator),
  SpaceObjectType: () => (/* reexport safe */ _types_types_js__WEBPACK_IMPORTED_MODULE_1__.SpaceObjectType),
  Star: () => (/* reexport safe */ _objects_index_js__WEBPACK_IMPORTED_MODULE_7__.Star),
  StateCovariance: () => (/* reexport safe */ _covariance_index_js__WEBPACK_IMPORTED_MODULE_17__.StateCovariance),
  StateVector: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.StateVector),
  Sun: () => (/* reexport safe */ _body_index_js__WEBPACK_IMPORTED_MODULE_8__.Sun),
  TAU: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.TAU),
  TEME: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.TEME),
  TimeStamped: () => (/* reexport safe */ _time_index_js__WEBPACK_IMPORTED_MODULE_3__.TimeStamped),
  Tle: () => (/* reexport safe */ _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__.Tle),
  Vector: () => (/* reexport safe */ _operations_operations_js__WEBPACK_IMPORTED_MODULE_6__.Vector),
  Vector3D: () => (/* reexport safe */ _operations_operations_js__WEBPACK_IMPORTED_MODULE_6__.Vector3D),
  ZoomValue: () => (/* reexport safe */ _types_types_js__WEBPACK_IMPORTED_MODULE_1__.ZoomValue),
  acoth: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.acoth),
  acsch: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.acsch),
  angularDiameter: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.angularDiameter),
  angularDistance: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.angularDistance),
  angularVelocityOfEarth: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.angularVelocityOfEarth),
  array2d: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.array2d),
  asec2rad: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.asec2rad),
  asech: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.asech),
  astronomicalUnit: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.astronomicalUnit),
  azel2uv: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.azel2uv),
  cKmPerMs: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.cKmPerMs),
  cKmPerSec: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.cKmPerSec),
  cMPerSec: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.cMPerSec),
  calcGmst: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.calcGmst),
  calcIncFromAz: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.calcIncFromAz),
  calcInertAz: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.calcInertAz),
  clamp: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.clamp),
  concat: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.concat),
  copySign: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.copySign),
  covariance: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.covariance),
  createCovarianceFromTle: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.createCovarianceFromTle),
  createSampleCovarianceFromTle: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.createSampleCovarianceFromTle),
  createVec: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.createVec),
  csch: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.csch),
  deg2rad: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.deg2rad),
  derivative: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.derivative),
  dopplerFactor: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.dopplerFactor),
  earthGravityParam: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.earthGravityParam),
  ecf2eci: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.ecf2eci),
  ecf2enu: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.ecf2enu),
  ecf2rae: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.ecf2rae),
  ecfRad2rae: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.ecfRad2rae),
  eci2ecf: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.eci2ecf),
  eci2lla: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.eci2lla),
  eci2rae: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.eci2rae),
  enu2rf: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.enu2rf),
  evalPoly: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.evalPoly),
  factorial: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.factorial),
  gamma: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.gamma),
  getDayOfYear: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.getDayOfYear),
  getDegLat: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.getDegLat),
  getDegLon: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.getDegLon),
  getRadLat: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.getRadLat),
  getRadLon: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.getRadLon),
  halfPi: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.halfPi),
  isLeapYear: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.isLeapYear),
  jacobian: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.jacobian),
  jday: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.jday),
  linearDistance: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.linearDistance),
  linearInterpolate: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.linearInterpolate),
  lla2ecef: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.lla2ecef),
  lla2ecf: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.lla2ecf),
  lla2eci: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.lla2eci),
  lla2sez: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.lla2sez),
  llaRad2ecf: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.llaRad2ecf),
  log10: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.log10),
  masec2rad: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.masec2rad),
  matchHalfPlane: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.matchHalfPlane),
  mean: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.mean),
  msec2sec: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.msec2sec),
  newtonM: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.newtonM),
  newtonNu: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.newtonNu),
  normalizeAngle: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.normalizeAngle),
  observationDerivative: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.observationDerivative),
  observationNoiseFromSigmas: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.observationNoiseFromSigmas),
  rad2deg: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rad2deg),
  radecToPosition: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.radecToPosition),
  radecToVelocity: () => (/* reexport safe */ _observation_index_js__WEBPACK_IMPORTED_MODULE_10__.radecToVelocity),
  rae2ecf: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2ecf),
  rae2eci: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2eci),
  rae2enu: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2enu),
  rae2raeOffBoresight: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2raeOffBoresight),
  rae2ruv: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2ruv),
  rae2sez: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.rae2sez),
  sec2day: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.sec2day),
  sec2deg: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.sec2deg),
  sec2min: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.sec2min),
  sech: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.sech),
  secondsPerDay: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.secondsPerDay),
  secondsPerSiderealDay: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.secondsPerSiderealDay),
  secondsPerWeek: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.secondsPerWeek),
  sez2rae: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.sez2rae),
  sign: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.sign),
  spaceObjType2Str: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.spaceObjType2Str),
  std: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.std),
  temp4: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.temp4),
  toPrecision: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.toPrecision),
  ttasec2rad: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.ttasec2rad),
  uv2azel: () => (/* reexport safe */ _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__.uv2azel),
  wrapAngle: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.wrapAngle),
  x2o3: () => (/* reexport safe */ _utils_index_js__WEBPACK_IMPORTED_MODULE_5__.x2o3)
});
/* ESM import */var _enums_index_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./enums/index.js */ "./src/engine/ootk/src/enums/index.ts");
/* ESM import */var _types_types_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./types/types.js */ "./src/engine/ootk/src/types/types.ts");
/* ESM import */var _interfaces_index_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./interfaces/index.js */ "./src/engine/ootk/src/interfaces/index.ts");
/* ESM import */var _time_index_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./time/index.js */ "./src/engine/ootk/src/time/index.ts");
/* ESM import */var _transforms_index_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./transforms/index.js */ "./src/engine/ootk/src/transforms/index.ts");
/* ESM import */var _utils_index_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./utils/index.js */ "./src/engine/ootk/src/utils/index.ts");
/* ESM import */var _operations_operations_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./operations/operations.js */ "./src/engine/ootk/src/operations/operations.ts");
/* ESM import */var _objects_index_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./objects/index.js */ "./src/engine/ootk/src/objects/index.ts");
/* ESM import */var _body_index_js__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./body/index.js */ "./src/engine/ootk/src/body/index.ts");
/* ESM import */var _coordinate_index_js__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./coordinate/index.js */ "./src/engine/ootk/src/coordinate/index.ts");
/* ESM import */var _observation_index_js__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./observation/index.js */ "./src/engine/ootk/src/observation/index.ts");
/* ESM import */var _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(/*! ./data/DataHandler.js */ "./src/engine/ootk/src/data/DataHandler.ts");
/* ESM import */var _sgp4_index_js__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(/*! ./sgp4/index.js */ "./src/engine/ootk/src/sgp4/index.ts");
/* ESM import */var _operations_index_js__WEBPACK_IMPORTED_MODULE_13__ = __webpack_require__(/*! ./operations/index.js */ "./src/engine/ootk/src/operations/index.ts");
/* ESM import */var _force_index_js__WEBPACK_IMPORTED_MODULE_14__ = __webpack_require__(/*! ./force/index.js */ "./src/engine/ootk/src/force/index.ts");
/* ESM import */var _propagator_index_js__WEBPACK_IMPORTED_MODULE_15__ = __webpack_require__(/*! ./propagator/index.js */ "./src/engine/ootk/src/propagator/index.ts");
/* ESM import */var _orbit_determination_index_js__WEBPACK_IMPORTED_MODULE_16__ = __webpack_require__(/*! ./orbit_determination/index.js */ "./src/engine/ootk/src/orbit_determination/index.ts");
/* ESM import */var _covariance_index_js__WEBPACK_IMPORTED_MODULE_17__ = __webpack_require__(/*! ./covariance/index.js */ "./src/engine/ootk/src/covariance/index.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






















}),
"./src/engine/ootk/src/objects/BaseObject.ts": 
/*!***************************************************!*\
  !*** ./src/engine/ootk/src/objects/BaseObject.ts ***!
  \***************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BaseObject: () => (BaseObject)
});
/* ESM import */var _types_types_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../types/types.js */ "./src/engine/ootk/src/types/types.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

class BaseObject {
    id;
    name;
    type;
    position; // Where is the object
    totalVelocity; // How fast is the object moving
    velocity; // How fast is the object moving
    active = true; // Is the object active
    constructor(info) {
        this.type = info.type ?? _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.UNKNOWN;
        this.name = info.name ?? 'Unknown';
        this.id = info.id ?? -1; // Default to -1 if no id is provided
        this.active = info.active ?? true;
        // Default to the center of the earth until position is calculated
        this.position = info.position ?? {
            x: 0,
            y: 0,
            z: 0,
        };
        // Default to 0 velocity until velocity is calculated
        this.velocity = info.velocity ?? {
            x: 0,
            y: 0,
            z: 0,
        };
        this.totalVelocity = Math.sqrt(this.velocity.x ** 2 + this.velocity.y ** 2 + this.velocity.z ** 2);
    }
    /**
     * Checks if the object is a satellite.
     * @returns True if the object is a satellite, false otherwise.
     */
    isSatellite() {
        return false;
    }
    /**
     * Checks if the object is a ground object.
     * @returns True if the object is a ground object, false otherwise.
     */
    isGroundObject() {
        return false;
    }
    /**
     * Returns whether the object is a sensor.
     * @returns True if the object is a sensor, false otherwise.
     */
    isSensor() {
        return false;
    }
    /**
     * Checks if the object is a marker.
     * @returns True if the object is a marker, false otherwise.
     */
    isMarker() {
        return false;
    }
    /**
     * Returns whether the object's position is static.
     * @returns True if the object is static, false otherwise.
     */
    isStatic() {
        return this.velocity.x === 0 && this.velocity.y === 0 && this.velocity.z === 0;
    }
    isPayload() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD;
    }
    isRocketBody() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.ROCKET_BODY;
    }
    isDebris() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.DEBRIS;
    }
    isStar() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.STAR;
    }
    isMissile() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.BALLISTIC_MISSILE;
    }
    isNotional() {
        return this.type === _types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.NOTIONAL;
    }
    getTypeString() {
        const typeToStringMap = {
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.UNKNOWN]: 'Unknown',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD]: 'Payload',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.ROCKET_BODY]: 'Rocket Body',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.DEBRIS]: 'Debris',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.SPECIAL]: 'Special',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.BALLISTIC_MISSILE]: 'Ballistic Missile',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.STAR]: 'Star',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.INTERGOVERNMENTAL_ORGANIZATION]: 'Intergovernmental Organization',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.SUBORBITAL_PAYLOAD_OPERATOR]: 'Suborbital Payload Operator',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD_OWNER]: 'Payload Owner',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.METEOROLOGICAL_ROCKET_LAUNCH_AGENCY_OR_MANUFACTURER]: 'Meteorological Rocket Launch Agency or Manufacturer',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD_MANUFACTURER]: 'Payload Manufacturer',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_AGENCY]: 'Launch Agency',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_SITE]: 'Launch Site',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_POSITION]: 'Launch Position',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_FACILITY]: 'Launch Facility',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.CONTROL_FACILITY]: 'Control Facility',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.GROUND_SENSOR_STATION]: 'Ground Sensor Station',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.OPTICAL]: 'Optical',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.MECHANICAL]: 'Mechanical',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PHASED_ARRAY_RADAR]: 'Phased Array Radar',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.OBSERVER]: 'Observer',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.BISTATIC_RADIO_TELESCOPE]: 'Bistatic Radio Telescope',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.COUNTRY]: 'Country',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_VEHICLE_MANUFACTURER]: 'Launch Vehicle Manufacturer',
            [_types_types_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.ENGINE_MANUFACTURER]: 'Engine Manufacturer',
        };
        return typeToStringMap[this.type] ?? 'Unknown';
    }
    /**
     * Validates a parameter value against a minimum and maximum value.
     * @param value - The value to be validated.
     * @param minValue - The minimum allowed value.
     * @param maxValue - The maximum allowed value.
     * @param errorMessage - The error message to be thrown if the value is invalid.
     */
    validateParameter(value, minValue, maxValue, errorMessage) {
        if (typeof minValue !== 'undefined' && minValue !== null && value < minValue) {
            throw new Error(errorMessage);
        }
        if (typeof maxValue !== 'undefined' && maxValue !== null && value > maxValue) {
            throw new Error(errorMessage);
        }
    }
}


}),
"./src/engine/ootk/src/objects/DetailedSatellite.ts": 
/*!**********************************************************!*\
  !*** ./src/engine/ootk/src/objects/DetailedSatellite.ts ***!
  \**********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DetailedSatellite: () => (DetailedSatellite)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _types_types_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../types/types.js */ "./src/engine/ootk/src/types/types.ts");
/* ESM import */var _Satellite_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Satellite.js */ "./src/engine/ootk/src/objects/Satellite.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



/**
 * Represents a detailed satellite object with launch, spacecraft, and operations details.
 */
class DetailedSatellite extends _Satellite_js__WEBPACK_IMPORTED_MODULE_2__.Satellite {
    configuration = '';
    country = '';
    dryMass = '';
    equipment = '';
    launchDate = '';
    launchMass = '';
    launchSite = '';
    launchPad = '';
    launchVehicle = '';
    lifetime = '';
    maneuver = '';
    manufacturer = '';
    mission = '';
    bus = '';
    motor = '';
    owner = '';
    payload = '';
    power = '';
    purpose = '';
    length = '';
    diameter = '';
    shape = '';
    span = '';
    user = '';
    source = '';
    vmag;
    rcs;
    altId = '';
    altName = '';
    status = _types_types_js__WEBPACK_IMPORTED_MODULE_1__.PayloadStatus.UNKNOWN;
    constructor(
    // TODO: Replace this intersection with a type alias
    info, options) {
        if (info.source === _main_js__WEBPACK_IMPORTED_MODULE_0__.CatalogSource.VIMPEL) {
            info = DetailedSatellite.setSccNumTo0_(info);
        }
        super(info, options);
        this.active ??= true;
        this.initSpaceCraftDetails_(info);
        this.length = info.length ?? '';
        this.diameter = info.diameter ?? '';
        this.source = info.source ?? '';
        this.vmag = info.vmag ?? null;
        this.rcs = info.rcs ?? null;
        this.altId = info.altId ?? '';
        this.altName = info.altName ?? '';
        this.initOperationDetails_(info);
        this.initLaunchDetails_(info);
        this.status = info.status ?? _types_types_js__WEBPACK_IMPORTED_MODULE_1__.PayloadStatus.UNKNOWN;
    }
    static setSccNumTo0_(info) {
        info.tle1 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle1, 2, '0');
        info.tle1 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle1, 3, '0');
        info.tle1 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle1, 4, '0');
        info.tle1 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle1, 5, '0');
        info.tle1 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle1, 6, '0');
        info.tle2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle2, 2, '0');
        info.tle2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle2, 3, '0');
        info.tle2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle2, 4, '0');
        info.tle2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle2, 5, '0');
        info.tle2 = _main_js__WEBPACK_IMPORTED_MODULE_0__.FormatTle.setCharAt(info.tle2, 6, '0');
        return info;
    }
    initSpaceCraftDetails_(info) {
        this.lifetime = info.lifetime ?? '';
        this.maneuver = info.maneuver ?? '';
        this.manufacturer = info.manufacturer ?? '';
        this.motor = info.motor ?? '';
        this.power = info.power ?? '';
        this.payload = info.payload ?? '';
        this.purpose = info.purpose ?? '';
        this.shape = info.shape ?? '';
        this.span = info.span ?? '';
        this.bus = info.bus ?? '';
        this.configuration = info.configuration ?? '';
        this.equipment = info.equipment ?? '';
        this.dryMass = info.dryMass ?? '';
    }
    initOperationDetails_(info) {
        this.mission = info.mission ?? '';
        this.user = info.user ?? '';
        this.owner = info.owner ?? '';
        this.country = info.country ?? '';
    }
    initLaunchDetails_(info) {
        this.launchDate = info.launchDate ?? '';
        this.launchMass = info.launchMass ?? '';
        this.launchSite = info.launchSite ?? '';
        this.launchPad = info.launchPad ?? '';
        this.launchVehicle = info.launchVehicle ?? '';
    }
    /**
     * Returns the launch details of the satellite.
     * @returns An object containing the launch date, launch mass, launch site, and launch vehicle of the satellite.
     */
    getLaunchDetails() {
        return {
            launchDate: this.launchDate,
            launchMass: this.launchMass,
            launchSite: this.launchSite,
            launchVehicle: this.launchVehicle,
        };
    }
    /**
     * Returns an object containing the details of the operations.
     * @returns An object containing the user, mission, owner, and country details.
     */
    getOperationsDetails() {
        return {
            user: this.user,
            mission: this.mission,
            owner: this.owner,
            country: this.country,
        };
    }
    /**
     * Returns the space craft details.
     * @returns The space craft details.
     */
    getSpaceCraftDetails() {
        return {
            lifetime: this.lifetime,
            maneuver: this.maneuver,
            manufacturer: this.manufacturer,
            motor: this.motor,
            power: this.power,
            payload: this.payload,
            purpose: this.purpose,
            shape: this.shape,
            span: this.span,
            configuration: this.configuration,
            equipment: this.equipment,
            dryMass: this.dryMass,
        };
    }
    clone() {
        return new DetailedSatellite(this);
    }
}


}),
"./src/engine/ootk/src/objects/DetailedSensor.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/objects/DetailedSensor.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DetailedSensor: () => (DetailedSensor)
});
/* ESM import */var _Sensor_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Sensor.js */ "./src/engine/ootk/src/objects/Sensor.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

class DetailedSensor extends _Sensor_js__WEBPACK_IMPORTED_MODULE_0__.Sensor {
    sensorId;
    objName;
    shortName;
    uiName;
    country;
    dwellTime;
    freqBand;
    commLinks;
    /** Is this sensor volumetric? */
    isVolumetric;
    /** The ideal zoom to see the sensor's full FOV */
    zoom;
    system;
    operator;
    url;
    constructor(info) {
        super(info);
        this.commLinks = info.commLinks ?? [];
        this.country = info.country;
        this.dwellTime = info.changeObjectInterval;
        this.freqBand = info.freqBand;
        this.isVolumetric = info.volume;
        this.objName = info.objName;
        this.operator = info.operator;
        this.sensorId = info.sensorId;
        this.shortName = info.shortName;
        this.system = info.system;
        this.uiName = info.uiName;
        this.url = info.url;
        this.zoom = info.zoom;
    }
    isStatic() {
        return true;
    }
}


}),
"./src/engine/ootk/src/objects/GroundObject.ts": 
/*!*****************************************************!*\
  !*** ./src/engine/ootk/src/objects/GroundObject.ts ***!
  \*****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  GroundObject: () => (GroundObject)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


class GroundObject extends _BaseObject_js__WEBPACK_IMPORTED_MODULE_1__.BaseObject {
    name = 'Unknown Ground Object';
    lat;
    lon;
    alt;
    constructor(info) {
        super(info);
        this.validateGroundObjectInputData_(info);
        this.name = info.name ?? this.name;
        this.lat = info.lat;
        this.lon = info.lon;
        this.alt = info.alt;
    }
    /**
     * Calculates the relative azimuth, elevation, and range between this GroundObject and a Satellite.
     * @param satellite The Satellite object.
     * @param date The date for which to calculate the RAE values. Defaults to the current date.
     * @returns The relative azimuth, elevation, and range values in kilometers and degrees.
     */
    rae(satellite, date = new Date()) {
        return satellite.rae(this, date);
    }
    /**
     * Calculates ECF position at a given time.
     * @variation optimized version of this.toGeodetic().toITRF().position;
     * @returns The ECF position vector of the ground object.
     */
    ecf() {
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.llaRad2ecf)(this.toGeodetic());
    }
    /**
     * Calculates the Earth-Centered Inertial (ECI) position vector of the ground object at a given date.
     * @variation optimzed version of this.toGeodetic().toITRF().toJ2000().position;
     * @param date The date for which to calculate the ECI position vector. Defaults to the current date.
     * @returns The ECI position vector of the ground object.
     */
    eci(date = new Date()) {
        const { gmst } = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.calcGmst)(date);
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.lla2eci)(this.toGeodetic(), gmst);
    }
    /**
     * Returns the latitude, longitude, and altitude of the GroundObject.
     * @returns The latitude, longitude, and altitude as an LlaVec3 object.
     */
    lla() {
        return {
            lat: this.lat,
            lon: this.lon,
            alt: this.alt,
        };
    }
    /**
     * Converts the latitude, longitude, and altitude of the GroundObject to radians and kilometers.
     * @variation optimized version of this.toGeodetic() without class instantiation for better performance and
     * serialization.
     * @returns An object containing the latitude, longitude, and altitude in
     * radians and kilometers.
     */
    llaRad() {
        return {
            lat: (this.lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
            lon: (this.lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
            alt: this.alt,
        };
    }
    get latRad() {
        return this.lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    }
    get lonRad() {
        return this.lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    }
    /**
     * Creates a GroundObject object from a Geodetic position.
     * @param geodetic The geodetic coordinates.
     * @returns A new GroundObject object.
     */
    static fromGeodetic(geodetic) {
        return new GroundObject({
            lat: geodetic.latDeg,
            lon: geodetic.lonDeg,
            alt: geodetic.alt,
        });
    }
    /**
     * Converts the ground position to geodetic coordinates.
     * @returns The geodetic coordinates.
     */
    toGeodetic() {
        return _main_js__WEBPACK_IMPORTED_MODULE_0__.Geodetic.fromDegrees(this.lat, this.lon, this.alt);
    }
    /**
     * Validates the input data for the GroundObject.
     * @param info - The GroundPositionParams object containing the latitude,
     * longitude, and altitude. @returns void
     */
    validateGroundObjectInputData_(info) {
        this.validateParameter(info.lat, -90, 90, 'Invalid latitude - must be between -90 and 90');
        this.validateParameter(info.lon, -180, 180, 'Invalid longitude - must be between -180 and 180');
        this.validateParameter(info.alt, 0, null, 'Invalid altitude - must be greater than 0');
    }
    isGroundObject() {
        switch (this.type) {
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.INTERGOVERNMENTAL_ORGANIZATION:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.SUBORBITAL_PAYLOAD_OPERATOR:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD_OWNER:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.METEOROLOGICAL_ROCKET_LAUNCH_AGENCY_OR_MANUFACTURER:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PAYLOAD_MANUFACTURER:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_VEHICLE_MANUFACTURER:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.ENGINE_MANUFACTURER:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_AGENCY:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_SITE:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.LAUNCH_POSITION:
                return true;
            default:
                return false;
        }
    }
}


}),
"./src/engine/ootk/src/objects/LandObject.ts": 
/*!***************************************************!*\
  !*** ./src/engine/ootk/src/objects/LandObject.ts ***!
  \***************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  LandObject: () => (LandObject)
});
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

class LandObject extends _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__.BaseObject {
    lat;
    lon;
    alt;
    country;
    Code;
    constructor(info) {
        super(info);
        this.lat = info.lat;
        this.lon = info.lon;
        this.alt = info.alt;
    }
    isLandObject() {
        return true;
    }
}


}),
"./src/engine/ootk/src/objects/Marker.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/objects/Marker.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Marker: () => (Marker)
});
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/* eslint-disable class-methods-use-this */
class Marker extends _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__.BaseObject {
    isMarker() {
        return true;
    }
}


}),
"./src/engine/ootk/src/objects/RfSensor.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/objects/RfSensor.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RfSensor: () => (RfSensor)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _DetailedSensor_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./DetailedSensor.js */ "./src/engine/ootk/src/objects/DetailedSensor.ts");
/**
 * @author Theodore Kruczek.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2022-2024 Theodore Kruczek Permission is
 * hereby granted, free of charge, to any person obtaining a copy of this
 * software and associated documentation files (the "Software"), to deal in the
 * Software without restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


class RfSensor extends _DetailedSensor_js__WEBPACK_IMPORTED_MODULE_1__.DetailedSensor {
    boresightAz;
    boresightEl;
    faces;
    beamwidth;
    constructor(info) {
        super(info);
        switch (info.type) {
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.BISTATIC_RADIO_TELESCOPE:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.MECHANICAL:
            case _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.PHASED_ARRAY_RADAR:
                break;
            default:
                throw new Error('Invalid sensor type');
        }
        this.boresightAz = info.boresightAz;
        this.boresightEl = info.boresightEl;
        if (info.boresightAz.length !== info.boresightEl.length) {
            throw new Error('Boresight azimuth and elevation arrays must be the same length');
        }
        this.faces = info.boresightAz.length;
        this.beamwidth = info.beamwidth;
    }
    /**
     * Converts azimuth and elevation angles to unit vector coordinates.
     * @param az - The azimuth angle in degrees.
     * @param el - The elevation angle in degrees.
     * @param face - The face number (optional).
     * @returns The unit vector coordinates.
     */
    uvFromAzEl(az, el, face) {
        const azRad = (az * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const elRad = (el * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
        const azDiff = (azRad - this.boresightAzRad(face ?? 0));
        const elDiff = (elRad - this.boresightElRad(face ?? 0));
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.azel2uv)(azDiff, elDiff, this.beamwidthRad);
    }
    /**
     * Converts the given UV coordinates to azimuth and elevation angles.
     * @param u - The U coordinate.
     * @param v - The V coordinate.
     * @param face - The face number for multi-faced sensors. (optional)
     * @returns An object containing the azimuth and elevation angles in degrees.
     * @throws Error if face number is not specified for multi-faced sensors.
     */
    azElFromUV(u, v, face) {
        if (!face && this.faces > 1) {
            throw new Error('Face number must be specified for multi-faced sensors');
        }
        const { az, el } = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.uv2azel)(u, v, this.beamwidthRad);
        return {
            az: ((az + this.boresightAz[face ?? 0]) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
            el: ((el + this.boresightEl[face ?? 0]) * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
        };
    }
    /**
     * Converts the boresight azimuth angle to radians.
     * @param face - The face number for multi-faced sensors. (Optional)
     * @returns The boresight azimuth angle in radians.
     * @throws An error if the face number is not specified for multi-faced sensors.
     */
    boresightAzRad(face) {
        if (!face && this.faces > 1) {
            throw new Error('Face number must be specified for multi-faced sensors');
        }
        return (this.boresightAz[face ?? 0] * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    }
    /**
     * Converts the boresight elevation angle of the sensor to radians.
     * @param face - The face number of the sensor (optional for single-faced sensors).
     * @returns The boresight elevation angle in radians.
     * @throws Error if the face number is not specified for multi-faced sensors.
     */
    boresightElRad(face) {
        if (!face && this.faces > 1) {
            throw new Error('Face number must be specified for multi-faced sensors');
        }
        return (this.boresightEl[face ?? 0] * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    }
    /**
     * Gets the beamwidth in radians.
     * @returns The beamwidth in radians.
     */
    get beamwidthRad() {
        return (this.beamwidth * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    }
}


}),
"./src/engine/ootk/src/objects/Satellite.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/objects/Satellite.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Satellite: () => (Satellite)
});
/* ESM import */var _coordinate_J2000_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../coordinate/J2000.js */ "./src/engine/ootk/src/coordinate/J2000.ts");
/* ESM import */var _coordinate_RIC_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../coordinate/RIC.js */ "./src/engine/ootk/src/coordinate/RIC.ts");
/* ESM import */var _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../coordinate/Tle.js */ "./src/engine/ootk/src/coordinate/Tle.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _observation_RAE_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../observation/RAE.js */ "./src/engine/ootk/src/observation/RAE.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ../time/EpochUTC.js */ "./src/engine/ootk/src/time/EpochUTC.ts");
/* ESM import */var _transforms_index_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ../transforms/index.js */ "./src/engine/ootk/src/transforms/index.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */











/**
 * Represents a satellite object with orbital information and methods for
 * calculating its position and other properties.
 */
class Satellite extends _BaseObject_js__WEBPACK_IMPORTED_MODULE_10__.BaseObject {
    apogee;
    argOfPerigee;
    bstar;
    eccentricity;
    epochDay;
    epochYear;
    inclination;
    intlDes;
    meanAnomaly;
    meanMoDev1;
    meanMoDev2;
    meanMotion;
    options;
    perigee;
    period;
    rightAscension;
    satrec;
    /** The satellite catalog number as listed in the TLE. */
    sccNum;
    /** The 5 digit alpha-numeric satellite catalog number. */
    sccNum5;
    /** The 6 digit numeric satellite catalog number. */
    sccNum6;
    tle1;
    tle2;
    /** The semi-major axis of the satellite's orbit. */
    semiMajorAxis;
    /** The semi-minor axis of the satellite's orbit. */
    semiMinorAxis;
    constructor(info, options) {
        super(info);
        if (info.tle1 && info.tle2) {
            this.parseTleAndUpdateOrbit_(info.tle1, info.tle2, info.sccNum);
        }
        else if (info.omm) {
            this.parseOmmAndUpdateOrbit_(info.omm);
        }
        else {
            throw new Error('tle1 and tle2 or omm must be provided to create a Satellite object.');
        }
        this.options = options ?? {
            notes: '',
        };
    }
    parseTleAndUpdateOrbit_(tle1, tle2, sccNum) {
        const tleData = _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.parse(tle1, tle2);
        this.tle1 = tle1;
        this.tle2 = tle2;
        this.sccNum = sccNum ?? tleData.satNum.toString();
        this.sccNum5 = _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.convert6DigitToA5(this.sccNum);
        this.sccNum6 = _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.convertA5to6Digit(this.sccNum5);
        this.intlDes = tleData.intlDes;
        this.epochYear = tleData.epochYear;
        this.epochDay = tleData.epochDay;
        this.meanMoDev1 = tleData.meanMoDev1;
        this.meanMoDev2 = tleData.meanMoDev2;
        this.bstar = tleData.bstar;
        this.inclination = tleData.inclination;
        this.rightAscension = tleData.rightAscension;
        this.eccentricity = tleData.eccentricity;
        this.argOfPerigee = tleData.argOfPerigee;
        this.meanAnomaly = tleData.meanAnomaly;
        this.meanMotion = tleData.meanMotion;
        this.period = tleData.period;
        this.semiMajorAxis = ((8681663.653 / this.meanMotion) ** (2 / 3));
        this.semiMinorAxis = (this.semiMajorAxis * Math.sqrt(1 - this.eccentricity ** 2));
        this.apogee = (this.semiMajorAxis * (1 + this.eccentricity) - 6371);
        this.perigee = (this.semiMajorAxis * (1 - this.eccentricity) - 6371);
        this.satrec = _main_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4.createSatrec(tle1, tle2);
    }
    parseOmmAndUpdateOrbit_(omm) {
        this.sccNum = omm.NORAD_CAT_ID.padStart(5, '0');
        this.sccNum5 = _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.convert6DigitToA5(omm.NORAD_CAT_ID);
        this.sccNum6 = _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.convertA5to6Digit(this.sccNum5);
        this.intlDes = omm.OBJECT_ID;
        const YYYY = omm.EPOCH.slice(0, 4);
        const MM = omm.EPOCH.slice(5, 7);
        const DD = omm.EPOCH.slice(8, 10);
        const hh = omm.EPOCH.slice(11, 13);
        const mm = omm.EPOCH.slice(14, 16);
        const ss = omm.EPOCH.slice(17, 23);
        const epochDateObj = Date.UTC(Number(YYYY), Number(MM) - 1, Number(DD), Number(hh), Number(mm), Number(ss));
        const dayOfYear = (epochDateObj - Date.UTC(Number(YYYY), 0, 0)) / 86400000;
        const ommParsed = {
            ...omm,
            epoch: {
                year: Number(YYYY),
                month: Number(MM),
                day: Number(DD),
                hour: Number(hh),
                minute: Number(mm),
                second: Number(ss),
                doy: dayOfYear,
            },
        };
        this.epochYear = parseInt(YYYY.slice(2, 4));
        this.epochDay = dayOfYear;
        this.meanMoDev1 = parseFloat(omm.MEAN_MOTION_DOT);
        this.meanMoDev2 = parseFloat(omm.MEAN_MOTION_DDOT);
        this.bstar = parseFloat(omm.BSTAR);
        this.inclination = parseFloat(omm.INCLINATION);
        this.rightAscension = parseFloat(omm.RA_OF_ASC_NODE);
        this.eccentricity = parseFloat(omm.ECCENTRICITY);
        this.argOfPerigee = parseFloat(omm.ARG_OF_PERICENTER);
        this.meanAnomaly = parseFloat(omm.MEAN_ANOMALY);
        this.meanMotion = parseFloat(omm.MEAN_MOTION);
        this.period = 1440 / this.meanMotion;
        this.semiMajorAxis = ((8681663.653 / this.meanMotion) ** (2 / 3));
        this.semiMinorAxis = (this.semiMajorAxis * Math.sqrt(1 - this.eccentricity ** 2));
        this.apogee = (this.semiMajorAxis * (1 + this.eccentricity) - 6371);
        this.perigee = (this.semiMajorAxis * (1 - this.eccentricity) - 6371);
        this.satrec = _main_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4.createSatrecFromOmm(ommParsed);
    }
    /**
     * Checks if the object is a satellite.
     * @returns True if the object is a satellite, false otherwise.
     */
    isSatellite() {
        return true;
    }
    /**
     * Returns whether the satellite is static or not.
     * @returns True if the satellite is static, false otherwise.
     */
    isStatic() {
        return false;
    }
    /**
     * Checks if the given SatelliteRecord object is valid by checking if its properties are all numbers.
     * @param satrec - The SatelliteRecord object to check.
     * @returns True if the SatelliteRecord object is valid, false otherwise.
     */
    static isValidSatrec(satrec) {
        if (isNaN(satrec.a) ||
            isNaN(satrec.am) ||
            isNaN(satrec.alta) ||
            isNaN(satrec.em) ||
            isNaN(satrec.mo) ||
            isNaN(satrec.ecco) ||
            isNaN(satrec.no)) {
            return false;
        }
        return true;
    }
    ageOfElset(nowInput, outputUnits = 'days') {
        return _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle.calcElsetAge(this.tle1, nowInput, outputUnits);
    }
    editTle(tle1, tle2, sccNum) {
        this.parseTleAndUpdateOrbit_(tle1, tle2, sccNum);
    }
    /**
     * Calculates the azimuth angle of the satellite relative to the given sensor at the specified date. If no date is
     * provided, the current time of the satellite is used.
     * @variation optimized
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the azimuth angle. Optional, defaults to the current date.
     * @returns The azimuth angle of the satellite relative to the given sensor at the specified date.
     */
    az(observer, date = new Date()) {
        const rae = this.rae(observer, date);
        if (!rae) {
            return null;
        }
        return (rae.az * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.RAD2DEG);
    }
    /**
     * Calculates the RAE (Range, Azimuth, Elevation) values for a given sensor and date. If no date is provided, the
     * current time is used.
     * @variation expanded
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the RAE values. Optional, defaults to the current date.
     * @returns The RAE values for the given sensor and date.
     */
    toRae(observer, date = new Date()) {
        const rae = this.rae(observer, date);
        if (!rae) {
            return null;
        }
        const rae2 = this.rae(observer, new Date(date.getTime() + 1000));
        if (!rae2) {
            return null;
        }
        const epoch = new _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_6__.EpochUTC(date.getTime() / 1000);
        const rangeRate = rae2.rng - rae.rng;
        const azimuthRate = rae2.az - rae.az;
        const elevationRate = rae2.el - rae.el;
        return new _observation_RAE_js__WEBPACK_IMPORTED_MODULE_4__.RAE(epoch, rae.rng, (rae.az * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.DEG2RAD), (rae.el * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.DEG2RAD), rangeRate, azimuthRate, elevationRate);
    }
    /**
     * Calculates ECF position at a given time.
     * @variation optimized
     * @param date - The date at which to calculate the ECF position. Optional, defaults to the current date.
     * @returns The ECF position at the specified date.
     */
    ecf(date = new Date()) {
        const { gmst } = Satellite.calculateTimeVariables(date);
        const eci = this.eci(date);
        if (!eci) {
            return null;
        }
        return (0,_transforms_index_js__WEBPACK_IMPORTED_MODULE_7__.eci2ecf)(eci.position, gmst);
    }
    /**
     * Calculates ECI position at a given time.
     * @variation optimized
     * @param date - The date at which to calculate the ECI position. Optional, defaults to the current date.
     * @param j - Julian date. Optional, defaults to null.
     * @param gmst - Greenwich Mean Sidereal Time. Optional, defaults to null.
     * @returns The ECI position at the specified date.
     */
    eci(date, j, gmst) {
        date ??= new Date();
        const { m } = Satellite.calculateTimeVariables(date, this.satrec, j, gmst);
        if (m === null) {
            return null;
        }
        const pv = _main_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4.propagate(this.satrec, m);
        if (!pv.position || !pv.velocity) {
            return null;
        }
        return pv;
    }
    /**
     * Calculates the J2000 coordinates for a given date. If no date is provided, the current time is used.
     * @variation expanded
     * @param date - The date for which to calculate the J2000 coordinates, defaults to the current date.
     * @returns The J2000 coordinates for the specified date.
     * @throws Error if propagation fails.
     */
    toJ2000(date = new Date()) {
        const { m } = Satellite.calculateTimeVariables(date, this.satrec);
        if (m === null) {
            throw new Error('Propagation failed!');
        }
        const pv = _main_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4.propagate(this.satrec, m);
        if (!pv.position) {
            throw new Error('Propagation failed!');
        }
        else {
            const p = pv.position;
            const v = pv.velocity;
            const epoch = new _time_EpochUTC_js__WEBPACK_IMPORTED_MODULE_6__.EpochUTC(date.getTime() / 1000);
            const pos = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_5__.Vector3D(p.x, p.y, p.z);
            const vel = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_5__.Vector3D(v.x, v.y, v.z);
            return new _coordinate_J2000_js__WEBPACK_IMPORTED_MODULE_0__.J2000(epoch, pos, vel);
        }
    }
    /**
     * Returns the elevation angle of the satellite as seen by the given sensor at the specified time.
     * @variation optimized
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the elevation angle. Optional, defaults to the current date.
     * @returns The elevation angle of the satellite as seen by the given sensor at the specified time.
     */
    el(observer, date = new Date()) {
        const rae = this.rae(observer, date);
        if (!rae) {
            return null;
        }
        return (rae.el * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.RAD2DEG);
    }
    /**
     * Calculates LLA position at a given time.
     * @variation optimized
     * @param date - The date at which to calculate the LLA position. Optional, defaults to the current date.
     * @param j - Julian date. Optional, defaults to null.
     * @param gmst - Greenwich Mean Sidereal Time. Optional, defaults to null.
     * @returns The LLA position at the specified date.
     */
    lla(date, j, gmst) {
        date ??= new Date();
        if (!j || !gmst) {
            const timeVar = Satellite.calculateTimeVariables(date, this.satrec);
            j = timeVar.j;
            gmst = timeVar.gmst;
        }
        const eci = this.eci(date, j, gmst);
        if (!eci) {
            return null;
        }
        const pos = eci.position;
        const lla = (0,_transforms_index_js__WEBPACK_IMPORTED_MODULE_7__.eci2lla)(pos, gmst);
        return lla;
    }
    /**
     * Converts the satellite's position to geodetic coordinates.
     * @variation expanded
     * @param date The date for which to calculate the geodetic coordinates. Defaults to the current date.
     * @returns The geodetic coordinates of the satellite.
     */
    toGeodetic(date = new Date()) {
        return this.toJ2000(date).toITRF().toGeodetic();
    }
    /**
     * Converts the satellite's position to the International Terrestrial Reference Frame (ITRF) at the specified date.
     * If no date is provided, the current date is used.
     * @variation expanded
     * @param date The date for which to convert the position. Defaults to the current date.
     * @returns The satellite's position in the ITRF at the specified date.
     */
    toITRF(date = new Date()) {
        return this.toJ2000(date).toITRF();
    }
    /**
     * Converts the current satellite's position to the Reference-Inertial-Celestial (RIC) frame
     * relative to the specified reference satellite at the given date.
     * @variation expanded
     * @param reference The reference satellite.
     * @param date The date for which to calculate the RIC frame. Defaults to the current date.
     * @returns The RIC frame representing the current satellite's position relative to the reference satellite.
     */
    toRIC(reference, date = new Date()) {
        return _coordinate_RIC_js__WEBPACK_IMPORTED_MODULE_1__.RIC.fromJ2000(this.toJ2000(date), reference.toJ2000(date));
    }
    /**
     * Converts the satellite object to a TLE (Two-Line Element) object.
     * @returns The TLE object representing the satellite.
     */
    toTle() {
        return new _coordinate_Tle_js__WEBPACK_IMPORTED_MODULE_2__.Tle(this.tle1, this.tle2);
    }
    /**
     * Converts the satellite's position to classical orbital elements.
     * @param date The date for which to calculate the classical elements. Defaults to the current date.
     * @returns The classical orbital elements of the satellite.
     */
    toClassicalElements(date = new Date()) {
        return this.toJ2000(date).toClassicalElements();
    }
    /**
     * Calculates the RAE (Range, Azimuth, Elevation) vector for a given sensor and time.
     * @variation optimized
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the RAE vector. Optional, defaults to the current date.
     * @param j - Julian date. Optional, defaults to null.
     * @param gmst - Greenwich Mean Sidereal Time. Optional, defaults to null.
     * @returns The RAE vector for the given sensor and time.
     */
    rae(observer, date, j, gmst) {
        date ??= new Date();
        gmst ??= Satellite.calculateTimeVariables(date, this.satrec).gmst;
        const eci = this.eci(date, j, gmst);
        if (!eci) {
            return null;
        }
        const ecf = (0,_transforms_index_js__WEBPACK_IMPORTED_MODULE_7__.eci2ecf)(eci.position, gmst);
        return (0,_transforms_index_js__WEBPACK_IMPORTED_MODULE_7__.ecf2rae)(observer, ecf);
    }
    /**
     * Returns the range of the satellite from the given sensor at the specified time.
     * @variation optimized
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the range. Optional, defaults to the current date.
     * @returns The range of the satellite from the given sensor at the specified time.
     */
    rng(observer, date = new Date()) {
        const rae = this.rae(observer, date);
        if (!rae) {
            return null;
        }
        return rae.rng;
    }
    /**
     * Applies the Doppler effect to the given frequency based on the observer's position and the date.
     * @param freq - The frequency to apply the Doppler effect to.
     * @param observer - The observer's position on the ground.
     * @param date - The date at which to calculate the Doppler effect. Optional, defaults to the current date.
     * @returns The frequency after applying the Doppler effect.
     */
    applyDoppler(freq, observer, date) {
        const doppler = this.dopplerFactor(observer, date);
        if (!doppler) {
            return null;
        }
        return freq * doppler;
    }
    /**
     * Calculates the Doppler factor for the satellite.
     * @param observer The observer's ground position.
     * @param date The optional date for which to calculate the Doppler factor. If not provided, the current date is used.
     * @returns The calculated Doppler factor.
     */
    dopplerFactor(observer, date) {
        const position = this.eci(date);
        if (!position) {
            return null;
        }
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_9__.dopplerFactor)(observer.eci(date), position.position, position.velocity);
    }
    /**
     * Calculates the time variables for a given date relative to the TLE epoch.
     * @param date Date to calculate
     * @param satrec Satellite orbital information
     * @param j Julian date
     * @param gmst Greenwich Mean Sidereal Time
     * @returns Time variables
     */
    static calculateTimeVariables(date, satrec, j, gmst) {
        j ??= (0,_transforms_index_js__WEBPACK_IMPORTED_MODULE_7__.jday)(date.getUTCFullYear(), date.getUTCMonth() + 1, date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds()) + date.getUTCMilliseconds() * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.MILLISECONDS_TO_DAYS;
        gmst ??= _main_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4.gstime(j);
        const m = satrec ? (j - satrec.jdsatepoch) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_8__.MINUTES_PER_DAY : null;
        return { gmst, m, j };
    }
}


}),
"./src/engine/ootk/src/objects/Sensor.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/objects/Sensor.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sensor: () => (Sensor)
});
/* ESM import */var _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/PassType.js */ "./src/engine/ootk/src/enums/PassType.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _types_types_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../types/types.js */ "./src/engine/ootk/src/types/types.ts");
/* ESM import */var _GroundObject_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./GroundObject.js */ "./src/engine/ootk/src/objects/GroundObject.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




class Sensor extends _GroundObject_js__WEBPACK_IMPORTED_MODULE_3__.GroundObject {
    minRng;
    minAz;
    minEl;
    maxRng;
    maxAz;
    maxEl;
    minRng2;
    minAz2;
    minEl2;
    maxRng2;
    maxAz2;
    maxEl2;
    constructor(info) {
        // If there is a sensor type verify it is valid
        if (info.type) {
            switch (info.type) {
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.OPTICAL:
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.MECHANICAL:
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.PHASED_ARRAY_RADAR:
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.OBSERVER:
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.BISTATIC_RADIO_TELESCOPE:
                case _types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.SHORT_TERM_FENCE:
                    break;
                default:
                    throw new Error('Invalid sensor type');
            }
        }
        super(info);
        this.validateSensorInputData_(info);
        this.minRng = info.minRng;
        this.minAz = info.minAz;
        this.minEl = info.minEl;
        this.maxRng = info.maxRng;
        this.maxAz = info.maxAz;
        this.maxEl = info.maxEl;
        this.minRng2 = info.minRng2;
        this.minAz2 = info.minAz2;
        this.minEl2 = info.minEl2;
        this.maxRng2 = info.maxRng2;
        this.maxAz2 = info.maxAz2;
        this.maxEl2 = info.maxEl2;
    }
    /**
     * Checks if the object is a sensor.
     * @returns True if the object is a sensor, false otherwise.
     */
    isSensor() {
        return true;
    }
    calculatePasses(planningInterval, sat, date = new Date()) {
        let isInViewLast = false;
        let maxElThisPass = 0;
        const msnPlanPasses = [];
        const startTime = date.getTime();
        for (let timeOffset = 0; timeOffset < planningInterval; timeOffset++) {
            const curTime = new Date(startTime + timeOffset * 1000);
            const rae = this.rae(sat, curTime);
            if (!rae) {
                continue;
            }
            const isInView = this.isRaeInFov(rae);
            if (timeOffset === 0) {
                // Propagate Backwards to get the previous pass
                const oldRae = this.rae(sat, new Date(date.getTime() - 1 * 1000));
                if (!oldRae) {
                    continue;
                }
                isInViewLast = this.isRaeInFov(oldRae);
            }
            const type = Sensor.getPassType_(isInView, isInViewLast);
            maxElThisPass = Math.max(maxElThisPass, rae.el);
            if (type === _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.ENTER || type === _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.EXIT) {
                const pass = {
                    type,
                    time: curTime,
                    az: rae.az,
                    el: rae.el,
                    rng: rae.rng,
                };
                // Only set maxEl for EXIT passes
                if (type === _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.EXIT) {
                    pass.maxElPass = maxElThisPass;
                }
                msnPlanPasses.push(pass);
                maxElThisPass = 0;
            }
            isInViewLast = isInView;
        }
        return msnPlanPasses;
    }
    /**
     * Checks if the given RAE vector is within the field of view of the sensor.
     * @param rae - The RAE vector to check.
     * @returns True if the RAE vector is within the field of view, false otherwise.
     */
    isRaeInFov(rae) {
        if (rae.el < this.minEl || rae.el > this.maxEl) {
            return false;
        }
        if (rae.rng < this.minRng || rae.rng > this.maxRng) {
            return false;
        }
        if (this.minAz > this.maxAz) {
            // North Facing Sensors
            if (rae.az < this.minAz && rae.az > this.maxAz) {
                return false;
            }
            // Normal Facing Sensors
        }
        else if (rae.az < this.minAz || rae.az > this.maxAz) {
            return false;
        }
        return true;
    }
    /**
     * Checks if a satellite is in the field of view (FOV) of the sensor.
     * @param sat - The satellite to check.
     * @param date - The date to use for the calculation. Defaults to the current date.
     * @returns A boolean indicating whether the satellite is in the FOV.
     */
    isSatInFov(sat, date = new Date()) {
        const rae = this.rae(sat, date);
        if (!rae) {
            return false;
        }
        return this.isRaeInFov(rae);
    }
    /**
     * Checks if the sensor is in deep space.
     * @returns True if the sensor is in deep space, false otherwise.
     */
    isDeepSpace() {
        return this.maxRng > 6000;
    }
    /**
     * Checks if the sensor is near Earth.
     * @returns True if the sensor is near Earth, false otherwise.
     */
    isNearEarth() {
        return this.maxRng <= 6000;
    }
    toJ2000(date = new Date()) {
        const gmst = (0,_main_js__WEBPACK_IMPORTED_MODULE_1__.calcGmst)(date).gmst;
        const position = (0,_main_js__WEBPACK_IMPORTED_MODULE_1__.lla2eci)(this.llaRad(), gmst);
        return new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(_main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC.fromDateTime(date), new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(position.x, position.y, position.z), new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(0, 0, 0));
    }
    /**
     * Returns the pass type based on the current and previous visibility states.
     * @param isInView - Indicates if the object is currently in view.
     * @param isInViewLast - Indicates if the object was in view in the previous state.
     * @returns The pass type.
     */
    static getPassType_(isInView, isInViewLast) {
        let type = _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.OUT_OF_VIEW;
        if (isInView && !isInViewLast) {
            type = _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.ENTER;
        }
        else if (!isInView && isInViewLast) {
            type = _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.EXIT;
        }
        else if (isInView && isInViewLast) {
            type = _enums_PassType_js__WEBPACK_IMPORTED_MODULE_0__.PassType.IN_VIEW;
        }
        return type;
    }
    /**
     * Validates the field of view (FOV) parameters of the sensor.
     * @param info - The sensor parameters.
     */
    validateFov_(info) {
        this.validateParameter(info.maxAz, 0, 360, 'Invalid maximum azimuth - must be between 0 and 360');
        this.validateParameter(info.minAz, 0, 360, 'Invalid maximum azimuth - must be between 0 and 360');
        this.validateParameter(info.maxEl, -15, 180, 'Invalid maximum elevation - must be between 0 and 180');
        this.validateParameter(info.minEl, -15, 90, 'Invalid minimum elevation - must be between 0 and 90');
        this.validateParameter(info.maxRng, 0, null, 'Invalid maximum range - must be greater than 0');
        this.validateParameter(info.minRng, 0, null, 'Invalid minimum range - must be greater than 0');
    }
    /**
     * Validates the field of view parameters for the sensor.
     * @param info - The sensor parameters.
     */
    validateFov2_(info) {
        this.validateParameter(info.maxAz2, 0, 360, 'Invalid maximum azimuth2 - must be between 0 and 360');
        this.validateParameter(info.minAz2, 0, 360, 'Invalid maximum azimuth2 - must be between 0 and 360');
        this.validateParameter(info.maxEl2, -15, 180, 'Invalid maximum elevation2 - must be between 0 and 180');
        this.validateParameter(info.minEl2, -15, 90, 'Invalid minimum elevation2 - must be between 0 and 90');
        this.validateParameter(info.maxRng2, 0, null, 'Invalid maximum range2 - must be greater than 0');
        this.validateParameter(info.minRng2, 0, null, 'Invalid minimum range2 - must be greater than 0');
    }
    /**
     * Validates the input data for the sensor.
     * @param info - The sensor parameters.
     */
    validateSensorInputData_(info) {
        this.validateLla_(info);
        this.validateFov_(info);
        if (info.minAz2 || info.maxAz2 || info.minEl2 || info.maxEl2 || info.minRng2 || info.maxRng2) {
            this.validateFov2_(info);
        }
    }
    /**
     * Validates the latitude, longitude, and altitude of a sensor.
     * @param info - The sensor parameters containing the latitude, longitude, and altitude.
     */
    validateLla_(info) {
        this.validateParameter(info.lat, -90, 90, 'Invalid latitude - must be between -90 and 90');
        this.validateParameter(info.lon, -180, 180, 'Invalid longitude - must be between -180 and 180');
        this.validateParameter(info.alt, 0, null, 'Invalid altitude - must be greater than 0');
    }
}


}),
"./src/engine/ootk/src/objects/Star.ts": 
/*!*********************************************!*\
  !*** ./src/engine/ootk/src/objects/Star.ts ***!
  \*********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Star: () => (Star)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


class Star extends _BaseObject_js__WEBPACK_IMPORTED_MODULE_1__.BaseObject {
    ra;
    dec;
    bf;
    h;
    pname;
    vmag;
    constructor(info) {
        super(info);
        this.type = _main_js__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.STAR;
        this.ra = info.ra;
        this.dec = info.dec;
        this.pname = info.pname ?? '';
        this.bf = info.bf ?? '';
        this.h = info.h ?? '';
        this.vmag = info.vmag;
    }
    eci(lla = { lat: 180, lon: 0, alt: 0 }, date = new Date()) {
        const rae = this.rae(lla, date);
        const { gmst } = Star.calculateTimeVariables_(date);
        // Arbitrary distance to enable using ECI coordinates
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.ecf2eci)((0,_main_js__WEBPACK_IMPORTED_MODULE_0__.rae2ecf)(rae, { lat: 0, lon: 0, alt: 0 }), gmst);
    }
    rae(lla = { lat: 180, lon: 0, alt: 0 }, date = new Date()) {
        const starPos = _main_js__WEBPACK_IMPORTED_MODULE_0__.Celestial.azEl(date, lla.lat, lla.lon, this.ra, this.dec);
        return { az: starPos.az, el: starPos.el, rng: 250000 };
    }
    static calculateTimeVariables_(date) {
        const j = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.jday)(date.getUTCFullYear(), date.getUTCMonth() + 1, date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds()) +
            date.getUTCMilliseconds() * _main_js__WEBPACK_IMPORTED_MODULE_0__.MILLISECONDS_TO_DAYS;
        const gmst = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4.gstime(j);
        return { gmst, j };
    }
}


}),
"./src/engine/ootk/src/objects/index.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/objects/index.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BaseObject: () => (/* reexport safe */ _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__.BaseObject),
  DetailedSatellite: () => (/* reexport safe */ _DetailedSatellite_js__WEBPACK_IMPORTED_MODULE_5__.DetailedSatellite),
  DetailedSensor: () => (/* reexport safe */ _DetailedSensor_js__WEBPACK_IMPORTED_MODULE_6__.DetailedSensor),
  GroundObject: () => (/* reexport safe */ _GroundObject_js__WEBPACK_IMPORTED_MODULE_1__.GroundObject),
  LandObject: () => (/* reexport safe */ _LandObject_js__WEBPACK_IMPORTED_MODULE_7__.LandObject),
  Marker: () => (/* reexport safe */ _Marker_js__WEBPACK_IMPORTED_MODULE_8__.Marker),
  RfSensor: () => (/* reexport safe */ _RfSensor_js__WEBPACK_IMPORTED_MODULE_9__.RfSensor),
  Satellite: () => (/* reexport safe */ _Satellite_js__WEBPACK_IMPORTED_MODULE_2__.Satellite),
  Sensor: () => (/* reexport safe */ _Sensor_js__WEBPACK_IMPORTED_MODULE_3__.Sensor),
  Star: () => (/* reexport safe */ _Star_js__WEBPACK_IMPORTED_MODULE_4__.Star)
});
/* ESM import */var _BaseObject_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./BaseObject.js */ "./src/engine/ootk/src/objects/BaseObject.ts");
/* ESM import */var _GroundObject_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./GroundObject.js */ "./src/engine/ootk/src/objects/GroundObject.ts");
/* ESM import */var _Satellite_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Satellite.js */ "./src/engine/ootk/src/objects/Satellite.ts");
/* ESM import */var _Sensor_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Sensor.js */ "./src/engine/ootk/src/objects/Sensor.ts");
/* ESM import */var _Star_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./Star.js */ "./src/engine/ootk/src/objects/Star.ts");
/* ESM import */var _DetailedSatellite_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./DetailedSatellite.js */ "./src/engine/ootk/src/objects/DetailedSatellite.ts");
/* ESM import */var _DetailedSensor_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./DetailedSensor.js */ "./src/engine/ootk/src/objects/DetailedSensor.ts");
/* ESM import */var _LandObject_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./LandObject.js */ "./src/engine/ootk/src/objects/LandObject.ts");
/* ESM import */var _Marker_js__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(/*! ./Marker.js */ "./src/engine/ootk/src/objects/Marker.ts");
/* ESM import */var _RfSensor_js__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(/*! ./RfSensor.js */ "./src/engine/ootk/src/objects/RfSensor.ts");












}),
"./src/engine/ootk/src/observation/ObservationUtils.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/observation/ObservationUtils.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  normalizeAngle: () => (normalizeAngle),
  observationDerivative: () => (observationDerivative),
  observationNoiseFromSigmas: () => (observationNoiseFromSigmas),
  radecToPosition: () => (radecToPosition),
  radecToVelocity: () => (radecToVelocity)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

const radecToPosition = (ra, dec, r) => {
    const ca = Math.cos(ra);
    const sa = Math.sin(ra);
    const cd = Math.cos(dec);
    const sd = Math.sin(dec);
    return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(r * cd * ca, r * cd * sa, r * sd);
};
const radecToVelocity = (ra, dec, r, raDot, decDot, rDot) => {
    const ca = Math.cos(ra);
    const sa = Math.sin(ra);
    const cd = Math.cos(dec);
    const sd = Math.sin(dec);
    return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(rDot * cd * ca - r * sd * ca * decDot - r * cd * sa * raDot, rDot * cd * sa - r * sd * sa * decDot + r * cd * ca * raDot, rDot * sd + r * cd * decDot);
};
const normalizeAngle = (a, b) => {
    const x = a - b;
    return Math.atan2(Math.sin(x), Math.cos(x));
};
const observationDerivative = (xh, xl, step, isAngle = false) => (isAngle ? normalizeAngle(xh, xl) : xh - xl) / step;
const observationNoiseFromSigmas = (sigmas) => {
    const n = sigmas.length;
    const result = Array.from({ length: n }, () => Array(n).fill(0.0));
    for (let i = 0; i < n; i++) {
        const s = sigmas[i];
        result[i][i] = 1 / (s * s);
    }
    return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix(result);
};


}),
"./src/engine/ootk/src/observation/PropagatorPairs.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/observation/PropagatorPairs.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  PropagatorPairs: () => (PropagatorPairs)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
class PropagatorPairs {
    posStep_;
    velStep_;
    constructor(posStep_, velStep_) {
        this.posStep_ = posStep_;
        this.velStep_ = velStep_;
        // Do nothing.
    }
    _high = Array(6).fill(null);
    _low = Array(6).fill(null);
    set(index, high, low) {
        this._high[index] = high;
        this._low[index] = low;
    }
    get(index) {
        return [this._high[index], this._low[index]];
    }
    /**
     * Get the step size at the provided index.
     * @param index The index.
     * @returns The step size.
     */
    step(index) {
        return index < 3 ? this.posStep_ : this.velStep_;
    }
}


}),
"./src/engine/ootk/src/observation/RAE.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/observation/RAE.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RAE: () => (RAE)
});
/* ESM import */var _coordinate_ITRF_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../coordinate/ITRF.js */ "./src/engine/ootk/src/coordinate/ITRF.ts");
/* ESM import */var _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../enums/AngularDistanceMethod.js */ "./src/engine/ootk/src/enums/AngularDistanceMethod.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable no-undefined */





// / Range, azimuth, and elevation.
class RAE {
    epoch;
    rng;
    azRad;
    elRad;
    rngRate;
    azRateRad;
    elRateRad;
    constructor(epoch, rng, azRad, elRad, 
    /** The range rate of the satellite relative to the observer in kilometers per second. */
    rngRate, 
    /** The azimuth rate of the satellite relative to the observer in radians per second. */
    azRateRad, 
    /** The elevation rate of the satellite relative to the observer in radians per second. */
    elRateRad) {
        this.epoch = epoch;
        this.rng = rng;
        this.azRad = azRad;
        this.elRad = elRad;
        this.rngRate = rngRate;
        this.azRateRad = azRateRad;
        this.elRateRad = elRateRad;
        // Do nothing
    }
    // / Create a new [Razel] object, using degrees for the angular values.
    static fromDegrees(epoch, range, azimuth, elevation, rangeRate, azimuthRate, elevationRate) {
        const azimuthRateRad = azimuthRate ? azimuthRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD : undefined;
        const elevationRateRad = elevationRate ? elevationRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD : undefined;
        return new RAE(epoch, range, (azimuth * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD), (elevation * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.DEG2RAD), rangeRate, azimuthRateRad, elevationRateRad);
    }
    /**
     * Create a [Razel] object from an inertial [state] and [site] vector.
     * @param state The inertial [state] vector.
     * @param site The observer [site] vector.
     * @returns A new [Razel] object.
     */
    static fromStateVector(state, site) {
        const stateEcef = state.toITRF();
        const siteEcef = site.toITRF();
        const po2 = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.halfPi;
        const r = stateEcef.position.subtract(siteEcef.position);
        const rDot = stateEcef.velocity;
        const geo = siteEcef.toGeodetic();
        const p = r.rotZ(geo.lon).rotY((po2 - geo.lat));
        const pDot = rDot.rotZ(geo.lon).rotY((po2 - geo.lat));
        const pS = p.x;
        const pE = p.y;
        const pZ = p.z;
        const pSDot = pDot.x;
        const pEDot = pDot.y;
        const pZDot = pDot.z;
        const pMag = p.magnitude();
        const pSEMag = Math.sqrt(pS * pS + pE * pE);
        const elevation = Math.asin(pZ / pMag);
        let azimuth;
        if (elevation !== po2) {
            azimuth = Math.atan2(-pE, pS) + Math.PI;
        }
        else {
            azimuth = Math.atan2(-pEDot, pSDot) + Math.PI;
        }
        const rangeRate = p.dot(pDot) / pMag;
        const azimuthRate = (pSDot * pE - pEDot * pS) / (pS * pS + pE * pE);
        const elevationRate = (pZDot - rangeRate * Math.sin(elevation)) / pSEMag;
        return new RAE(state.epoch, pMag, (azimuth % _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.TAU), elevation, rangeRate, azimuthRate, elevationRate);
    }
    /**
     * Gets the azimuth in degrees.
     * @returns The azimuth in degrees.
     */
    get az() {
        return this.azRad * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG;
    }
    /**
     * Gets the elevation angle in degrees.
     * @returns The elevation angle in degrees.
     */
    get el() {
        return this.elRad * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG;
    }
    /**
     * Gets the azimuth rate in degrees per second.
     * @returns The azimuth rate in degrees per second, or undefined if it is not available.
     */
    get azRate() {
        return this.azRateRad ? this.azRateRad * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG : undefined;
    }
    /**
     * Gets the elevation rate in degrees per second.
     * @returns The elevation rate in degrees per second, or undefined if the elevation rate is not set.
     */
    get elRate() {
        return this.elRateRad ? this.elRateRad * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG : undefined;
    }
    toString() {
        return [
            '[RazEl]',
            `  Epoch:     ${this.epoch}`,
            `  Azimuth:   ${this.az.toFixed(4)}°`,
            `  Elevation: ${this.el.toFixed(4)}°`,
            `  Range:     ${this.rng.toFixed(3)} km`,
        ].join('\n');
    }
    /**
     * Return the position relative to the observer [site].
     *
     * An optional azimuth [az] _(rad)_ and elevation [el] _(rad)_ value can be
     * passed to override the values contained in this observation.
     * @param site The observer [site].
     * @param azRad Azimuth _(rad)_.
     * @param elRad Elevation _(rad)_.
     * @returns A [Vector3D] object.
     */
    position(site, azRad, elRad) {
        const ecef = site.toITRF();
        const geo = ecef.toGeodetic();
        const po2 = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.halfPi;
        const newAz = azRad ?? this.azRad;
        const newEl = elRad ?? this.elRad;
        const sAz = Math.sin(newAz);
        const cAz = Math.cos(newAz);
        const sEl = Math.sin(newEl);
        const cEl = Math.cos(newEl);
        const pSez = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D((-this.rng * cEl * cAz), (this.rng * cEl * sAz), (this.rng * sEl));
        const rEcef = pSez
            .rotY(-(po2 - geo.lat))
            .rotZ(-geo.lon)
            .add(ecef.position);
        return new _coordinate_ITRF_js__WEBPACK_IMPORTED_MODULE_0__.ITRF(this.epoch, rEcef, _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D.origin).toJ2000().position;
    }
    /**
     * Convert this observation into a [J2000] state vector.
     *
     * This will throw an error if the [rangeRate], [elevationRate], or
     * [azimuthRate] are not defined.
     * @param site The observer [site].
     * @returns A [J2000] state vector.
     */
    toStateVector(site) {
        // If the rates are not defined then assume stationary
        this.rngRate ??= 0;
        this.elRateRad ??= 0;
        this.azRateRad ??= 0;
        const ecef = site.toITRF();
        const geo = ecef.toGeodetic();
        const po2 = _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.halfPi;
        const sAz = Math.sin(this.azRad);
        const cAz = Math.cos(this.azRad);
        const sEl = Math.sin(this.elRad);
        const cEl = Math.cos(this.elRad);
        const pSez = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D((-this.rng * cEl * cAz), (this.rng * cEl * sAz), (this.rng * sEl));
        const pDotSez = new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_2__.Vector3D((-this.rngRate * cEl * cAz +
            this.rng * sEl * cAz * this.elRateRad +
            this.rng * cEl * sAz * this.azRateRad), (this.rngRate * cEl * sAz -
            this.rng * sEl * sAz * this.elRateRad +
            this.rng * cEl * cAz * this.azRateRad), (this.rngRate * sEl + this.rng * cEl * this.elRateRad));
        const pEcef = pSez.rotY(-(po2 - geo.lat)).rotZ(-geo.lon);
        const pDotEcef = pDotSez
            .rotY(-(po2 - geo.lat))
            .rotZ(-geo.lon);
        const rEcef = pEcef.add(ecef.position);
        return new _coordinate_ITRF_js__WEBPACK_IMPORTED_MODULE_0__.ITRF(this.epoch, rEcef, pDotEcef).toJ2000();
    }
    /**
     * Calculate the angular distance _(rad)_ between this and another [Razel]
     * object.
     * @param razel The other [Razel] object.
     * @param method The angular distance method to use.
     * @returns The angular distance _(rad)_.
     */
    angle(razel, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Cosine) {
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_4__.angularDistance)(this.azRad, this.elRad, razel.azRad, razel.elRad, method);
    }
    /**
     * Calculate the angular distance _(°)_ between this and another [Razel]
     * object.
     * @param razel The other [Razel] object.
     * @param method The angular distance method to use.
     * @returns The angular distance _(°)_.
     */
    angleDegrees(razel, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Cosine) {
        return this.angle(razel, method) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_3__.RAD2DEG;
    }
}


}),
"./src/engine/ootk/src/observation/RadecGeocentric.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/observation/RadecGeocentric.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RadecGeocentric: () => (RadecGeocentric)
});
/* ESM import */var _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/AngularDistanceMethod.js */ "./src/engine/ootk/src/enums/AngularDistanceMethod.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./ObservationUtils.js */ "./src/engine/ootk/src/observation/ObservationUtils.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




/**
 * Represents a geocentric right ascension and declination observation.
 *
 * In geocentric coordinates, observations are considered from the Earth's center. This approach simplifies calculations
 * for distant celestial objects, as it assumes a uniform observation point that ignores the observer's specific
 * location on Earth.
 */
class RadecGeocentric {
    epoch;
    rightAscension;
    declination;
    range;
    rightAscensionRate;
    declinationRate;
    rangeRate;
    constructor(epoch, rightAscension, declination, range, rightAscensionRate, declinationRate, rangeRate) {
        this.epoch = epoch;
        this.rightAscension = rightAscension;
        this.declination = declination;
        this.range = range;
        this.rightAscensionRate = rightAscensionRate;
        this.declinationRate = declinationRate;
        this.rangeRate = rangeRate;
        // Nothing to do here.
    }
    /**
     * Creates a RadecGeocentric object from the given parameters in degrees.
     * @param epoch - The epoch in UTC.
     * @param rightAscensionDegrees - The right ascension in degrees.
     * @param declinationDegrees - The declination in degrees.
     * @param range - The range in kilometers (optional).
     * @param rightAscensionRateDegrees - The right ascension rate in degrees per second (optional).
     * @param declinationRateDegrees - The declination rate in degrees per second (optional).
     * @param rangeRate - The range rate in kilometers per second (optional).
     * @returns A new RadecGeocentric object.
     */
    static fromDegrees(epoch, rightAscensionDegrees, declinationDegrees, range, rightAscensionRateDegrees, declinationRateDegrees, rangeRate) {
        const rightAscensionRate = rightAscensionRateDegrees
            ? rightAscensionRateDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.DEG2RAD
            : null;
        const declinationRate = declinationRateDegrees ? declinationRateDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.DEG2RAD : null;
        return new RadecGeocentric(epoch, rightAscensionDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.DEG2RAD, declinationDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.DEG2RAD, range, rightAscensionRate, declinationRate, rangeRate);
    }
    /**
     * Creates a RadecGeocentric object from a state vector in J2000 coordinates.
     * @param state - The J2000 state vector.
     * @returns A new RadecGeocentric object.
     */
    static fromStateVector(state) {
        const rI = state.position.x;
        const rJ = state.position.y;
        const rK = state.position.z;
        const vI = state.velocity.x;
        const vJ = state.velocity.y;
        const vK = state.velocity.z;
        const rMag = state.position.magnitude();
        const declination = Math.asin(rK / rMag);
        const rIJMag = Math.sqrt(rI * rI + rJ * rJ);
        let rightAscension;
        if (rIJMag !== 0) {
            rightAscension = Math.atan2(rJ, rI);
        }
        else {
            rightAscension = Math.atan2(vJ, vI);
        }
        const rangeRate = state.position.dot(state.velocity) / rMag;
        const rightAscensionRate = (vI * rJ - vJ * rI) / (-(rJ * rJ) - rI * rI);
        const declinationRate = (vK - rangeRate * (rK / rMag)) / rIJMag;
        return new RadecGeocentric(state.epoch, rightAscension % _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.TAU, declination, rMag, rightAscensionRate, declinationRate, rangeRate);
    }
    /**
     * Gets the right ascension in degrees.
     * @returns The right ascension in degrees.
     */
    get rightAscensionDegrees() {
        return this.rightAscension * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG;
    }
    /**
     * Gets the declination in degrees.
     * @returns The declination in degrees.
     */
    get declinationDegrees() {
        return this.declination * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG;
    }
    /**
     * Gets the right ascension rate in degrees per second.
     * @returns The right ascension rate in degrees per second, or null if it is not available.
     */
    get rightAscensionRateDegrees() {
        return this.rightAscensionRate ? this.rightAscensionRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG : null;
    }
    /**
     * Gets the rate of change of declination in degrees per second.
     * @returns The rate of change of declination in degrees per second, or null if not available.
     */
    get declinationRateDegrees() {
        return this.declinationRate ? this.declinationRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG : null;
    }
    /**
     * Calculates the position vector in geocentric coordinates.
     * @param range - The range in kilometers (optional). If not provided, it uses the default range or 1.0 kilometer.
     * @returns The position vector in geocentric coordinates.
     */
    position(range) {
        const r = range ?? this.range ?? 1.0;
        return (0,_ObservationUtils_js__WEBPACK_IMPORTED_MODULE_3__.radecToPosition)(this.rightAscension, this.declination, r);
    }
    /**
     * Calculates the velocity vector of the celestial object.
     * @param range - The range of the celestial object in kilometers. If not provided, it uses the stored range value.
     * @param rangeRate - The range rate of the celestial object in kilometers per second.
     * If not provided, it uses the stored range rate value.
     * @returns The velocity vector of the celestial object in kilometers per second.
     * @throws Error if the right ascension rate or declination rate is missing.
     */
    velocity(range, rangeRate) {
        if (!this.rightAscensionRate || !this.declinationRate) {
            throw new Error('Velocity unsolvable, missing ra/dec rates.');
        }
        const r = range ?? this.range ?? 1.0;
        const rd = rangeRate ?? this.rangeRate ?? 0.0;
        return (0,_ObservationUtils_js__WEBPACK_IMPORTED_MODULE_3__.radecToVelocity)(this.rightAscension, this.declination, r, this.rightAscensionRate, this.declinationRate, rd);
    }
    /**
     * Calculates the angular distance between two celestial coordinates (RA and Dec).
     * @param radec - The celestial coordinates to compare with.
     * @param method - The method to use for calculating the angular distance. Default is `AngularDistanceMethod.Cosine`.
     * @returns The angular distance between the two celestial coordinates in radians.
     */
    angle(radec, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDistanceMethod.Cosine) {
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_2__.angularDistance)(this.rightAscension, this.declination, radec.rightAscension, radec.declination, method);
    }
    /**
     * Calculates the angle in degrees between two RadecGeocentric objects.
     * @param radec - The RadecGeocentric object to calculate the angle with.
     * @param method - The method to use for calculating the angular distance. Default is AngularDistanceMethod.Cosine.
     * @returns The angle in degrees.
     */
    angleDegrees(radec, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDistanceMethod.Cosine) {
        return this.angle(radec, method) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG;
    }
}


}),
"./src/engine/ootk/src/observation/RadecTopocentric.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/observation/RadecTopocentric.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RadecTopocentric: () => (RadecTopocentric)
});
/* ESM import */var _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/AngularDistanceMethod.js */ "./src/engine/ootk/src/enums/AngularDistanceMethod.ts");
/* ESM import */var _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../operations/Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./ObservationUtils.js */ "./src/engine/ootk/src/observation/ObservationUtils.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */





/**
 * Represents a topocentric right ascension and declination observation.
 *
 * Topocentric coordinates take into account the observer's exact location on the Earth's surface. This model is crucial
 * for precise measurements of local astronomical events and nearby celestial objects, where the observer's latitude,
 * longitude, and altitude can significantly affect the observed position due to parallax. Topocentric coordinates are
 * particularly important for observations of the Moon, planets, and artificial satellites.
 */
class RadecTopocentric {
    epoch;
    rightAscension;
    declination;
    range;
    rightAscensionRate;
    declinationRate;
    rangeRate;
    constructor(epoch, rightAscension, declination, range, rightAscensionRate, declinationRate, rangeRate) {
        this.epoch = epoch;
        this.rightAscension = rightAscension;
        this.declination = declination;
        this.range = range;
        this.rightAscensionRate = rightAscensionRate;
        this.declinationRate = declinationRate;
        this.rangeRate = rangeRate;
        // Nothing to do here.
    }
    /**
     * Create a new RadecTopocentric object, using degrees for the angular values.
     * @param epoch UTC epoch.
     * @param rightAscensionDegrees Right-ascension in degrees.
     * @param declinationDegrees Declination in degrees.
     * @param range Range in km.
     * @param rightAscensionRateDegrees Right-ascension rate in degrees per second.
     * @param declinationRateDegrees Declination rate in degrees per second.
     * @param rangeRate Range rate in km/s.
     * @returns A new RadecTopocentric object.
     */
    static fromDegrees(epoch, rightAscensionDegrees, declinationDegrees, range, rightAscensionRateDegrees, declinationRateDegrees, rangeRate) {
        const rightAscensionRate = rightAscensionRateDegrees
            ? rightAscensionRateDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD
            : null;
        const declinationRate = declinationRateDegrees
            ? declinationRateDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD
            : null;
        return new RadecTopocentric(epoch, rightAscensionDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD, declinationDegrees * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD, range, rightAscensionRate, declinationRate, rangeRate);
    }
    /**
     * Create a new RadecTopocentric object from a J2000 state vector.
     * @param state Inertial state vector.
     * @param site Site vector.
     * @returns A new RadecTopocentric object.
     */
    static fromStateVector(state, site) {
        const p = state.position.subtract(site.position);
        const pI = p.x;
        const pJ = p.y;
        const pK = p.z;
        const pMag = p.magnitude();
        const declination = Math.asin(pK / pMag);
        const pDot = state.velocity.subtract(site.velocity);
        const pIDot = pDot.x;
        const pJDot = pDot.y;
        const pKDot = pDot.z;
        const pIJMag = Math.sqrt(pI * pI + pJ * pJ);
        let rightAscension;
        if (pIJMag !== 0) {
            rightAscension = Math.atan2(pJ, pI);
        }
        else {
            rightAscension = Math.atan2(pJDot, pIDot);
        }
        const rangeRate = p.dot(pDot) / pMag;
        const rightAscensionRate = (pIDot * pJ - pJDot * pI) / (-(pJ * pJ) - pI * pI);
        const declinationRate = (pKDot - rangeRate * Math.sin(declination)) / pIJMag;
        return new RadecTopocentric(state.epoch, rightAscension % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU, declination, pMag, rightAscensionRate, declinationRate, rangeRate);
    }
    /**
     * Gets the right ascension in degrees.
     * @returns The right ascension in degrees.
     */
    get rightAscensionDegrees() {
        return this.rightAscension * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.RAD2DEG;
    }
    /**
     * Gets the declination in degrees.
     * @returns The declination in degrees.
     */
    get declinationDegrees() {
        return this.declination * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.RAD2DEG;
    }
    /**
     * Gets the right ascension rate in degrees per second.
     * @returns The right ascension rate in degrees per second, or null if it is not available.
     */
    get rightAscensionRateDegrees() {
        return this.rightAscensionRate
            ? this.rightAscensionRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.RAD2DEG
            : null;
    }
    /**
     * Gets the rate of change of declination in degrees per second.
     * @returns The rate of change of declination in degrees per second, or null if the declination rate is not defined.
     */
    get declinationRateDegrees() {
        return this.declinationRate
            ? this.declinationRate * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.RAD2DEG
            : null;
    }
    /**
     * Return the position relative to the observer site.
     *
     * An optional range value can be passed to override the value contained in this observation.
     * @param site Observer site.
     * @param range Range in km.
     * @returns A Vector3D object.
     */
    position(site, range) {
        const r = range ?? this.range ?? 1.0;
        return (0,_ObservationUtils_js__WEBPACK_IMPORTED_MODULE_4__.radecToPosition)(this.rightAscension, this.declination, r).add(site.position);
    }
    /**
     * Return the velocity relative to the observer site.
     *
     * An optional range and rangeRate value can be passed to override the values contained in this observation.
     * @param site Observer site.
     * @param range Range in km.
     * @param rangeRate Range rate in km/s.
     * @returns A Vector3D object.
     */
    velocity(site, range, rangeRate) {
        if (!this.rightAscensionRate || !this.declinationRate) {
            throw new Error('Velocity unsolvable, missing ra/dec rates.');
        }
        const r = range ?? this.range ?? 1.0;
        const rd = rangeRate ?? this.rangeRate ?? 0.0;
        return (0,_ObservationUtils_js__WEBPACK_IMPORTED_MODULE_4__.radecToVelocity)(this.rightAscension, this.declination, r, this.rightAscensionRate, this.declinationRate, rd).add(site.velocity);
    }
    /**
     * Calculates the line of sight vector in the topocentric coordinate system.
     * The line of sight vector points from the observer's location towards the celestial object.
     * @returns The line of sight vector as a Vector3D object.
     */
    lineOfSight() {
        const ca = Math.cos(this.rightAscension);
        const cd = Math.cos(this.declination);
        const sa = Math.sin(this.rightAscension);
        const sd = Math.sin(this.declination);
        return new _operations_Vector3D_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(cd * ca, cd * sa, sd);
    }
    /**
     * Calculate the angular distance between this and another RadecTopocentric object.
     * @param radec - The other RadecTopocentric object.
     * @param method - The angular distance method to use.
     * @returns The angular distance.
     */
    angle(radec, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDistanceMethod.Cosine) {
        return (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_3__.angularDistance)(this.rightAscension, this.declination, radec.rightAscension, radec.declination, method);
    }
    /**
     * Calculate the angular distance between this and another RadecTopocentric object.
     * @param radec - The other RadecTopocentric object.
     * @param method - The angular distance method to use.
     * @returns The angular distance
     */
    angleDegrees(radec, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDistanceMethod.Cosine) {
        return this.angle(radec, method) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.RAD2DEG;
    }
}


}),
"./src/engine/ootk/src/observation/index.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/observation/index.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RAE: () => (/* reexport safe */ _RAE_js__WEBPACK_IMPORTED_MODULE_3__.RAE),
  RadecGeocentric: () => (/* reexport safe */ _RadecGeocentric_js__WEBPACK_IMPORTED_MODULE_1__.RadecGeocentric),
  RadecTopocentric: () => (/* reexport safe */ _RadecTopocentric_js__WEBPACK_IMPORTED_MODULE_2__.RadecTopocentric),
  normalizeAngle: () => (/* reexport safe */ _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__.normalizeAngle),
  observationDerivative: () => (/* reexport safe */ _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__.observationDerivative),
  observationNoiseFromSigmas: () => (/* reexport safe */ _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__.observationNoiseFromSigmas),
  radecToPosition: () => (/* reexport safe */ _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__.radecToPosition),
  radecToVelocity: () => (/* reexport safe */ _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__.radecToVelocity)
});
/* ESM import */var _ObservationUtils_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./ObservationUtils.js */ "./src/engine/ootk/src/observation/ObservationUtils.ts");
/* ESM import */var _RadecGeocentric_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./RadecGeocentric.js */ "./src/engine/ootk/src/observation/RadecGeocentric.ts");
/* ESM import */var _RadecTopocentric_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./RadecTopocentric.js */ "./src/engine/ootk/src/observation/RadecTopocentric.ts");
/* ESM import */var _RAE_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./RAE.js */ "./src/engine/ootk/src/observation/RAE.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






}),
"./src/engine/ootk/src/operations/BoxMuller.ts": 
/*!*****************************************************!*\
  !*** ./src/engine/ootk/src/operations/BoxMuller.ts ***!
  \*****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BoxMuller: () => (BoxMuller)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Box-Muller random Gaussian number generator.
class BoxMuller {
    _index = 0;
    _cache = new Float64Array(2);
    // / Mean value.
    mu;
    // / Standard deviation.
    sigma;
    // / Uniform random number generator.
    rand;
    /**
     * Create a new [BoxMuller] object with mean [mu], standard deviation
     * [sigma], and [seed] number.
     * @param mu Mean value.
     * @param sigma Standard deviation.
     * @param seed Random seed.
     */
    constructor(mu, sigma, seed = 0) {
        this.mu = mu;
        this.sigma = sigma;
        this.rand = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Random(seed);
        this._generate();
    }
    // / Refill the cache with random Gaussian numbers.
    _generate() {
        this._index = 0;
        const u1 = this.rand.nextFloat();
        const u2 = this.rand.nextFloat();
        const mag = this.sigma * Math.sqrt(-2.0 * Math.log(u1));
        this._cache[0] = mag * Math.cos(_main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * u2) + this.mu;
        this._cache[1] = mag * Math.sin(_main_js__WEBPACK_IMPORTED_MODULE_0__.TAU * u2) + this.mu;
    }
    /**
     * Generate a gaussian number, with mean [mu] and standard
     * deviation [sigma].
     * @returns A gaussian number.
     */
    nextGauss() {
        if (this._index > 1) {
            this._generate();
        }
        const result = this._cache[this._index];
        this._index++;
        return result;
    }
    /**
     * Generate a [Vector] of gaussian numbers, with mean [mu] and standard
     * deviation [sigma].
     * @param n Number of gaussian numbers to generate.
     * @returns A [Vector] of gaussian numbers.
     */
    gaussVector(n) {
        const result = new Float64Array(n);
        for (let i = 0; i < n; i++) {
            result[i] = this.nextGauss();
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(result);
    }
}


}),
"./src/engine/ootk/src/operations/EulerAngles.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/operations/EulerAngles.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EulerAngles: () => (EulerAngles)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Class containing Euler angles.
class EulerAngles {
    // / Roll component _(rad)_.
    roll;
    // / Pitch component _(rad)_.
    pitch;
    // / Yaw component _(rad)_.
    yaw;
    /**
     * Create a new EulerAngles object from roll, pitch, and yaw angles in radians.
     * @param roll Roll angle in radians.
     * @param pitch Pitch angle in radians.
     * @param yaw Yaw angle in radians.
     */
    constructor(roll, pitch, yaw) {
        this.roll = roll;
        this.pitch = pitch;
        this.yaw = yaw;
    }
    /**
     * Create a new EulerAngles object from roll, pitch, and yaw angles.
     * @param rollDeg Roll angle in degrees.
     * @param pitchDeg Pitch angle in degrees.
     * @param yawDeg Yaw angle in degrees.
     * @returns EulerAngles object.
     */
    static fromDegrees(rollDeg, pitchDeg, yawDeg) {
        const roll = rollDeg * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
        const pitch = pitchDeg * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
        const yaw = yawDeg * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
        return new EulerAngles(roll, pitch, yaw);
    }
    /**
     * Create a new EulerAngles object from 3-2-1 ordered direction cosine matrix c.
     * @param c 3-2-1 ordered direction cosine matrix.
     * @returns EulerAngles object.
     */
    static fromDcm321(c) {
        const roll = Math.atan(c.elements[1][2] / c.elements[2][2]);
        const pitch = -Math.asin(c.elements[0][2]);
        const yaw = Math.atan(c.elements[0][1] / c.elements[0][0]);
        return new EulerAngles(roll, pitch, yaw);
    }
    /**
     * Gets the roll angle in degrees.
     * @returns The roll angle in degrees.
     */
    get rollDegrees() {
        return this.roll * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Gets the pitch angle in degrees.
     * @returns The pitch angle in degrees.
     */
    get pitchDegrees() {
        return this.pitch * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Gets the yaw angle in degrees.
     * @returns The yaw angle in degrees.
     */
    get yawDegrees() {
        return this.yaw * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Gets the roll angle in radians.
     * @returns The roll angle in radians.
     */
    get phi() {
        return this.roll;
    }
    /**
     * Gets the pitch angle in radians.
     * @returns The pitch angle in radians.
     */
    get theta() {
        return this.pitch;
    }
    /**
     * Gets the yaw angle in radians.
     * @returns The yaw angle in radians.
     */
    get psi() {
        return this.yaw;
    }
    /**
     * Gets the roll component in degrees.
     * @returns The roll component in degrees.
     */
    get phiDegrees() {
        return this.phi * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Gets the pitch component in degrees.
     * @returns The pitch component in degrees.
     */
    get thetaDegrees() {
        return this.theta * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Gets the yaw component in degrees.
     * @returns The yaw component in degrees.
     */
    get psiDegrees() {
        return this.psi * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    /**
     * Returns a string representation of the Euler angles.
     * @param precision The number of decimal places to include in the string representation. Default is 6.
     * @returns A string representation of the Euler angles.
     */
    toString(precision = 6) {
        const rollStr = this.rollDegrees.toFixed(precision);
        const pitchStr = this.pitchDegrees.toFixed(precision);
        const yawStr = this.yawDegrees.toFixed(precision);
        return `Euler(roll: ${rollStr}°, pitch: ${pitchStr}°, yaw: ${yawStr}°)`;
    }
    /**
     * Calculates the Direction Cosine Matrix (DCM) using the 3-2-1 Euler angles convention.
     * @returns The calculated DCM as a Matrix object.
     */
    dcm321() {
        const sPhi = Math.sin(this.phi);
        const cPhi = Math.cos(this.phi);
        const sTheta = Math.sin(this.theta);
        const cTheta = Math.cos(this.theta);
        const sPsi = Math.sin(this.psi);
        const cPsi = Math.cos(this.psi);
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [cTheta * cPsi, cTheta * sPsi, -sTheta],
            [sPhi * sTheta * cPsi - cPhi * sPsi, sPhi * sTheta * sPsi + cPhi * cPsi, sPhi * cTheta],
            [cPhi * sTheta * cPsi + sPhi * sPsi, cPhi * sTheta * sPsi - sPhi * cPsi, cPhi * cTheta],
        ]);
    }
    /**
     * Rotates a 3D vector using a 3-2-1 Euler angle sequence.
     * @param v The vector to rotate.
     * @returns The rotated vector.
     */
    rotateVector321(v) {
        return this.dcm321().multiplyVector3D(v);
    }
}


}),
"./src/engine/ootk/src/operations/Matrix.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/operations/Matrix.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Matrix: () => (Matrix)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * A matrix is a rectangular array of numbers or other mathematical objects for
 * which operations such as addition and multiplication are defined.
 */
class Matrix {
    elements;
    rows;
    columns;
    constructor(elements) {
        this.elements = elements;
        this.rows = elements.length;
        this.columns = (elements[0]).length;
    }
    /**
     * Creates a matrix with all elements set to zero.
     * @param rows - The number of rows in the matrix.
     * @param columns - The number of columns in the matrix.
     * @returns A matrix with all elements set to zero.
     */
    static allZeros(rows, columns) {
        return this.fill(rows, columns, 0.0);
    }
    /**
     * Creates a new Matrix with the specified number of rows and columns, filled
     * with the specified value.
     * @param rows The number of rows in the matrix.
     * @param columns The number of columns in the matrix.
     * @param value The value to fill the matrix with. Default is 0.0.
     * @returns A new Matrix filled with the specified value.
     */
    static fill(rows, columns, value = 0.0) {
        const elements = [];
        for (let i = 0; i < rows; i++) {
            elements[i] = [];
            for (let j = 0; j < columns; j++) {
                (elements[i])[j] = value;
            }
        }
        return new Matrix(elements);
    }
    /**
     * Creates a rotation matrix around the X-axis.
     * @param theta - The angle of rotation in radians.
     * @returns The rotation matrix.
     */
    static rotX(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const result = Matrix.zero(3, 3);
        (result.elements[0])[0] = 1.0;
        (result.elements[1])[1] = cosT;
        (result.elements[1])[2] = sinT;
        (result.elements[2])[1] = -sinT;
        (result.elements[2])[2] = cosT;
        return result;
    }
    /**
     * Creates a rotation matrix around the y-axis.
     * @param theta - The angle of rotation in radians.
     * @returns The rotation matrix.
     */
    static rotY(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const result = Matrix.zero(3, 3);
        (result.elements[0])[0] = cosT;
        (result.elements[0])[2] = -sinT;
        (result.elements[1])[1] = 1.0;
        (result.elements[2])[0] = sinT;
        (result.elements[2])[2] = cosT;
        return result;
    }
    /**
     * Creates a rotation matrix around the Z-axis.
     * @param theta The angle of rotation in radians.
     * @returns The rotation matrix.
     */
    static rotZ(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const result = Matrix.zero(3, 3);
        (result.elements[0])[0] = cosT;
        (result.elements[0])[1] = sinT;
        (result.elements[1])[0] = -sinT;
        (result.elements[1])[1] = cosT;
        (result.elements[2])[2] = 1.0;
        return result;
    }
    /**
     * Creates a zero matrix with the specified number of rows and columns.
     * @param rows The number of rows in the matrix.
     * @param columns The number of columns in the matrix.
     * @returns A new Matrix object representing the zero matrix.
     */
    static zero(rows, columns) {
        const elements = [];
        for (let i = 0; i < rows; i++) {
            elements[i] = [];
            for (let j = 0; j < columns; j++) {
                (elements[i])[j] = 0.0;
            }
        }
        return new Matrix(elements);
    }
    /**
     * Creates an identity matrix of the specified dimension.
     * @param dimension The dimension of the identity matrix.
     * @returns The identity matrix.
     */
    static identity(dimension) {
        const elements = [];
        for (let i = 0; i < dimension; i++) {
            elements[i] = [];
            for (let j = 0; j < dimension; j++) {
                (elements[i])[j] = i === j ? 1.0 : 0.0;
            }
        }
        return new Matrix(elements);
    }
    /**
     * Creates a diagonal matrix with the given diagonal elements.
     * @param d - An array of diagonal elements.
     * @returns A new Matrix object representing the diagonal matrix.
     */
    static diagonal(d) {
        const dimension = d.length;
        const elements = [];
        for (let i = 0; i < dimension; i++) {
            elements[i] = [];
            for (let j = 0; j < dimension; j++) {
                (elements[i])[j] = i === j ? d[i] : 0.0;
            }
        }
        return new Matrix(elements);
    }
    /**
     * Adds the elements of another matrix to this matrix and returns the result.
     * @param m - The matrix to be added.
     * @returns The resulting matrix after addition.
     */
    add(m) {
        const result = Matrix.zero(this.rows, this.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                (result.elements[i] ?? [])[j] = (this.elements[i]?.[j]) + (m.elements[i]?.[j]);
            }
        }
        return result;
    }
    /**
     * Subtracts the elements of another matrix from this matrix.
     * @param m - The matrix to subtract.
     * @returns A new matrix containing the result of the subtraction.
     */
    subtract(m) {
        const result = Matrix.zero(this.rows, this.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                (result.elements[i] ?? [])[j] = (this.elements[i]?.[j]) - (m.elements[i]?.[j]);
            }
        }
        return result;
    }
    /**
     * Scales the matrix by multiplying each element by a scalar value.
     * @param n - The scalar value to multiply each element by.
     * @returns A new Matrix object representing the scaled matrix.
     */
    scale(n) {
        const result = Matrix.zero(this.rows, this.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                (result.elements[i] ?? [])[j] = (this.elements[i]?.[j]) * n;
            }
        }
        return result;
    }
    /**
     * Negates the matrix by scaling it by -1.
     * @returns The negated matrix.
     */
    negate() {
        return this.scale(-1);
    }
    /**
     * Multiplies this matrix with another matrix.
     * @param m The matrix to multiply with.
     * @returns The resulting matrix.
     */
    multiply(m) {
        const result = Matrix.zero(this.rows, m.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < m.columns; j++) {
                for (let k = 0; k < this.columns; k++) {
                    ((result.elements[i])[j]) +=
                        (this.elements[i]?.[k]) * (m.elements[k]?.[j]);
                }
            }
        }
        return result;
    }
    /**
     * Computes the outer product of this matrix with another matrix.
     * @param m - The matrix to compute the outer product with.
     * @returns The resulting matrix.
     */
    outerProduct(m) {
        const result = Matrix.zero(this.rows, this.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                (result.elements[i])[j] = (this.elements[i]?.[j]) * (m.elements[i]?.[j]);
            }
        }
        return result;
    }
    /**
     * Multiplies the matrix by a vector.
     * @param v The vector to multiply by.
     * @returns A new vector representing the result of the multiplication.
     */
    multiplyVector(v) {
        const result = [];
        for (let i = 0; i < this.rows; i++) {
            let total = 0.0;
            for (let j = 0; j < this.columns; j++) {
                total += (this.elements[i]?.[j]) * (v.elements[j]);
            }
            result[i] = total;
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(result);
    }
    /**
     * Multiplies a 3D vector by the matrix.
     * @template T - The type of the vector elements.
     * @param v - The 3D vector to multiply.
     * @returns The resulting 3D vector after multiplication.
     */
    multiplyVector3D(v) {
        const result = [];
        for (let i = 0; i < this.rows; i++) {
            let total = 0.0;
            for (let j = 0; j < this.columns; j++) {
                switch (j) {
                    case 0:
                        total += (this.elements[i]?.[j]) * v.x;
                        break;
                    case 1:
                        total += (this.elements[i]?.[j]) * v.y;
                        break;
                    case 2:
                        total += (this.elements[i]?.[j]) * v.z;
                        break;
                    default:
                        break;
                }
            }
            result[i] = total;
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D((result[0]), (result[1]), (result[2]));
    }
    /**
     * Returns a new Matrix object where each element is the reciprocal of the
     * corresponding element in the current matrix. If an element in the current
     * matrix is zero, the corresponding element in the output matrix will also be
     * zero.
     * @returns A new Matrix object representing the reciprocal of the current
     * matrix.
     */
    reciprocal() {
        const output = Matrix.zero(this.rows, this.columns);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                if ((this.elements[i]?.[j]) !== 0) {
                    (output.elements[i])[j] = 1 / (this.elements[i]?.[j]);
                }
            }
        }
        return output;
    }
    /**
     * Transposes the matrix by swapping rows with columns.
     * @returns A new Matrix object representing the transposed matrix.
     */
    transpose() {
        const result = Matrix.zero(this.columns, this.rows);
        for (let i = 0; i < this.rows; i++) {
            for (let j = 0; j < this.columns; j++) {
                (result.elements[j])[i] = (this.elements[i])[j];
            }
        }
        return result;
    }
    /**
     * Performs the Cholesky decomposition on the matrix.
     * @returns A new Matrix object representing the Cholesky decomposition of the
     * original matrix.
     */
    cholesky() {
        const result = Matrix.zero(this.rows, this.rows);
        for (let i = 0; i < this.rows; i++) {
            for (let k = 0; k < i + 1; k++) {
                let total = 0.0;
                for (let j = 0; j < k; j++) {
                    total += (result.elements[i]?.[j]) * (result.elements[k]?.[j]);
                }
                (result.elements[i])[k] =
                    i === k
                        ? Math.sqrt((this.elements[i]?.[i]) - total)
                        : (1 / (result.elements[k]?.[k])) * ((this.elements[i]?.[k]) - total);
            }
        }
        return result;
    }
    /**
     * Swaps two rows in the matrix.
     * @param i - The index of the first row.
     * @param j - The index of the second row.
     */
    _swapRows(i, j) {
        if (i === j) {
            return;
        }
        const tmp = this.elements[i];
        this.elements[i] = this.elements[j];
        this.elements[j] = tmp;
    }
    /**
     * Converts the matrix to reduced row echelon form using the Gaussian
     * elimination method. This method modifies the matrix in-place.
     */
    toReducedRowEchelonForm_() {
        for (let lead = 0, row = 0; row < this.rows && lead < this.columns; ++row, ++lead) {
            let i = row;
            while ((this.elements[i]?.[lead]) === 0) {
                if (++i === this.rows) {
                    i = row;
                    if (++lead === this.columns) {
                        return;
                    }
                }
            }
            this._swapRows(i, row);
            if ((this.elements[row]?.[lead]) !== 0) {
                const f = this.elements[row]?.[lead];
                for (let column = 0; column < this.columns; ++column) {
                    (this.elements[row])[column] /= f;
                }
            }
            for (let j = 0; j < this.rows; ++j) {
                if (j === row) {
                    continue;
                }
                const f = (this.elements[j]?.[lead]);
                for (let column = 0; column < this.columns; ++column) {
                    ((this.elements[j])[column]) -= f * (this.elements[row]?.[column]);
                }
            }
        }
    }
    /**
     * Calculates the inverse of the matrix.
     * @returns The inverse of the matrix.
     */
    inverse() {
        const tmp = Matrix.zero(this.rows, this.columns * 2);
        for (let row = 0; row < this.rows; ++row) {
            for (let column = 0; column < this.columns; ++column) {
                (tmp.elements[row])[column] = (this.elements[row])[column];
            }
            (tmp.elements[row])[row + this.columns] = 1.0;
        }
        tmp.toReducedRowEchelonForm_();
        const inv = Matrix.zero(this.rows, this.columns);
        for (let row = 0; row < this.rows; ++row) {
            for (let column = 0; column < this.columns; ++column) {
                ((inv.elements[row])[column]) = (tmp.elements[row]?.[column + this.columns]);
            }
        }
        return inv;
    }
}


}),
"./src/engine/ootk/src/operations/Quaternion.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/operations/Quaternion.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Quaternion: () => (Quaternion)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

class Quaternion {
    x;
    y;
    z;
    w;
    constructor(x, y, z, w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
    static zero = new Quaternion(0, 0, 0, 0);
    static one = new Quaternion(0, 0, 0, 1);
    static xAxis = new Quaternion(1, 0, 0, 0);
    static yAxis = new Quaternion(0, 1, 0, 0);
    static zAxis = new Quaternion(0, 0, 1, 0);
    toString(precision = 8) {
        const xStr = this.x.toFixed(precision);
        const yStr = this.y.toFixed(precision);
        const zStr = this.z.toFixed(precision);
        const wStr = this.w.toFixed(precision);
        return `Q(x: ${xStr}, y: ${yStr}, z: ${zStr}, w: ${wStr})`;
    }
    positivePolar() {
        return this.w >= 0 ? this.normalize() : this.negate().normalize();
    }
    magnitudeSquared() {
        return this.w * this.w + this.x * this.x + this.y * this.y + this.z * this.z;
    }
    magnitude() {
        return Math.sqrt(this.magnitudeSquared());
    }
    scale(n) {
        return new Quaternion(n * this.x, n * this.y, n * this.z, n * this.w);
    }
    negate() {
        return this.scale(-1);
    }
    normalize() {
        const m = this.magnitude();
        if (m === 0) {
            return Quaternion.zero;
        }
        return this.scale(1 / m);
    }
    conjugate() {
        return new Quaternion(-this.x, -this.y, -this.z, this.w);
    }
    inverse() {
        return this.conjugate().scale(1 / this.magnitudeSquared());
    }
    add(q) {
        return new Quaternion(this.x + q.x, this.y + q.y, this.z + q.z, this.w + q.w);
    }
    subtract(q) {
        return new Quaternion(this.x - q.x, this.y - q.y, this.z - q.z, this.w - q.w);
    }
    addReal(n) {
        return new Quaternion(this.x, this.y, this.z, this.w + n);
    }
    multiply(q) {
        const mx = this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y;
        const my = this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x;
        const mz = this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w;
        const mw = this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z;
        return new Quaternion(mx, my, mz, mw);
    }
    dot(q) {
        return this.x * q.x + this.y * q.y + this.z * q.z + this.w * q.w;
    }
    rotateVector(v) {
        const q = this.multiply(new Quaternion(v.x, v.y, v.z, 0)).multiply(this.conjugate());
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector([q.x, q.y, q.z]);
    }
    rotateVector3D(v) {
        const q = this.multiply(new Quaternion(v.x, v.y, v.z, 0)).multiply(this.conjugate());
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(q.x, q.y, q.z);
    }
    lerp(q, t) {
        const f = 1.0 - t;
        return new Quaternion(f * this.x + t * q.x, f * this.y + t * q.y, f * this.z + t * q.z, f * this.w + t * q.w).positivePolar();
    }
    slerp(q, t) {
        let qp = q;
        let dotP = this.dot(qp);
        if (dotP < 0) {
            dotP = -dotP;
            qp = qp.negate();
        }
        if (dotP > 0.9995) {
            return this.lerp(qp, t);
        }
        const theta = Math.acos(dotP);
        const sinTheta = Math.sin(theta);
        const f1 = Math.sin((1.0 - t) * theta) / sinTheta;
        const f2 = Math.sin(t * theta) / sinTheta;
        return new Quaternion(f1 * this.x + f2 * qp.x, f1 * this.y + f2 * qp.y, f1 * this.z + f2 * qp.z, f1 * this.w + f2 * qp.w).positivePolar();
    }
    toVector3D() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(this.x, this.y, this.z);
    }
    angle(q) {
        const c = this.multiply(q.conjugate()).normalize();
        return 2 * Math.atan2(c.toVector3D().magnitude(), c.w);
    }
    geodesicAngle(q) {
        const p = this.dot(q);
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.wrapAngle)(Math.acos(2 * p * p - 1.0));
    }
    distance(q) {
        const m01 = this.subtract(q).magnitude();
        const p01 = this.add(q).magnitude();
        return m01 < p01 ? m01 : p01;
    }
    delta(qTo) {
        return this.inverse().multiply(qTo);
    }
    toDirectionCosineMatrix() {
        const w2 = this.w * this.w;
        const x2 = this.x * this.x;
        const y2 = this.y * this.y;
        const z2 = this.z * this.z;
        const m = [
            [w2 + x2 - y2 - z2, 2.0 * (this.x * this.y + this.z * this.w), 2.0 * (this.x * this.z - this.y * this.w)],
            [2.0 * (this.x * this.y - this.z * this.w), w2 - x2 + y2 - z2, 2.0 * (this.y * this.z + this.x * this.w)],
            [2.0 * (this.x * this.z + this.y * this.w), 2.0 * (this.y * this.z - this.x * this.w), w2 - x2 - y2 + z2],
        ];
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix(m);
    }
    toRotationMatrix() {
        return this.toDirectionCosineMatrix().transpose();
    }
    vectorAngle(observer, target, forward) {
        const delta = target.subtract(observer);
        const transform = this.toDirectionCosineMatrix().multiplyVector3D(delta);
        return forward.angle(transform);
    }
    kinematics(angularVelocity) {
        const wPrime = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector([0, angularVelocity.x, angularVelocity.y, angularVelocity.z]);
        const qMat = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [this.x, this.w, -this.z, this.y],
            [this.y, this.z, this.w, -this.x],
            [this.z, -this.y, this.x, this.w],
            [this.w, -this.x, -this.y, -this.z],
        ]);
        const result = qMat.multiplyVector(wPrime).scale(0.5).elements;
        return new Quaternion(result[0], result[1], result[2], result[3]);
    }
}


}),
"./src/engine/ootk/src/operations/Random.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/operations/Random.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Random: () => (Random)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/**
 * A generator of random bool, int, or double values.
 *
 * The default implementation supplies a stream of pseudo-random bits that are not suitable for cryptographic purposes.
 */
class Random {
    _seed;
    constructor(seed = 0) {
        this._seed = seed;
    }
    nextFloat(max = 1) {
        this._seed = (this._seed * 9301 + 49297) % 233280;
        return (this._seed / 233280) * max;
    }
    /**
     * To create a non-negative random integer uniformly distributed in the range from 0,
     * inclusive, to max, exclusive, use nextInt(int max).
     * @param max The bound on the random number to be returned. Must be positive.
     * @returns A pseudorandom, uniformly distributed int value between 0 (inclusive) and the specified value (exclusive).
     */
    nextInt(max = 1) {
        return Math.round(this.nextFloat(max) * max);
    }
    nextBool() {
        return this.nextFloat() > 0.5;
    }
}


}),
"./src/engine/ootk/src/operations/RandomGaussianSource.ts": 
/*!****************************************************************!*\
  !*** ./src/engine/ootk/src/operations/RandomGaussianSource.ts ***!
  \****************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RandomGaussianSource: () => (RandomGaussianSource)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _BoxMuller_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./BoxMuller.js */ "./src/engine/ootk/src/operations/BoxMuller.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


class RandomGaussianSource {
    boxMuller_;
    constructor(seed = 0) {
        this.boxMuller_ = new _BoxMuller_js__WEBPACK_IMPORTED_MODULE_1__.BoxMuller(0, 1, seed);
    }
    nextGauss() {
        return this.boxMuller_.nextGauss();
    }
    gaussVector(n) {
        if (n < 1) {
            throw new Error('n must be greater than 0');
        }
        const result = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector([this.nextGauss()]);
        for (let i = 0; i < n; i++) {
            if (i > 0) {
                result.add(new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector([this.nextGauss()]));
            }
        }
        return result;
    }
    gaussSphere(radius = 1.0) {
        return this.gaussVector(3).toVector3D(0).normalize().scale(radius);
    }
}


}),
"./src/engine/ootk/src/operations/Vector.ts": 
/*!**************************************************!*\
  !*** ./src/engine/ootk/src/operations/Vector.ts ***!
  \**************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Vector: () => (Vector)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * A Vector is a mathematical object that has both magnitude and direction.
 */
class Vector {
    elements;
    /**
     * The length of the vector.
     */
    length;
    /**
     * Represents a 3-dimensional vector.
     */
    static origin3 = new Vector([0, 0, 0]);
    /**
     * Represents a vector with all elements set to zero.
     */
    static origin6 = new Vector([0, 0, 0, 0, 0, 0]);
    /**
     * Represents the x-axis vector.
     */
    static xAxis = new Vector([1, 0, 0]);
    /**
     * Represents the y-axis vector.
     */
    static yAxis = new Vector([0, 1, 0]);
    /**
     * Represents the z-axis vector.
     */
    static zAxis = new Vector([0, 0, 1]);
    /**
     * Represents a vector pointing along the negative x-axis.
     */
    static xAxisNeg = new Vector([-1, 0, 0]);
    /**
     * Represents a vector pointing along the negative y-axis.
     */
    static yAxisNeg = new Vector([0, -1, 0]);
    /**
     * Represents a vector pointing along the negative z-axis.
     */
    static zAxisNeg = new Vector([0, 0, -1]);
    constructor(elements) {
        this.elements = elements;
        this.length = elements.length;
    }
    /**
     * Creates a zero vector of the specified length.
     * @param length The length of the vector.
     * @returns A new Vector object representing the zero vector.
     */
    static zero(length) {
        return new Vector(new Array(length).fill(0));
    }
    /**
     * Creates a new Vector with the specified length, filled with the specified
     * value.
     * @param length The length of the new Vector.
     * @param value The value to fill the Vector with.
     * @returns A new Vector filled with the specified value.
     */
    static filled(length, value) {
        return new Vector(new Array(length).fill(value));
    }
    /**
     * Creates a new Vector instance from an array of elements.
     * @param elements - The array of elements to create the Vector from.
     * @returns A new Vector instance.
     */
    static fromList(elements) {
        return new Vector(elements);
    }
    /**
     * Returns a string representation of the vector.
     * @param fixed - The number of digits to appear after the decimal point.
     * Defaults to -1.
     * @returns A string representation of the vector.
     */
    toString(fixed = -1) {
        if (fixed < 0) {
            return `[${this.elements.join(', ')}]`;
        }
        const output = this.elements.map((e) => e.toFixed(fixed));
        return `[${output.join(', ')}]`;
    }
    /**
     * Returns a string representation of the x value of the vector.
     * @returns A string representation of the x value of the vector.
     */
    get x() {
        return this.elements[0];
    }
    /**
     * Returns a string representation of the y value of the vector.
     * @returns A string representation of the y value of the vector.
     */
    get y() {
        return this.elements[1];
    }
    /**
     * Returns a string representation of the z value of the vector.
     * @returns A string representation of the z value of the vector.
     */
    get z() {
        return this.elements[2];
    }
    /**
     * Converts the vector elements to an array.
     * @returns An array containing the vector elements.
     */
    toList() {
        return Array.from(this.elements);
    }
    /**
     * Converts the vector to a Float64Array.
     * @returns The vector as a Float64Array.
     */
    toArray() {
        return new Float64Array(this.elements);
    }
    /**
     * Calculates the magnitude of the vector.
     * @returns The magnitude of the vector.
     */
    magnitude() {
        let total = 0;
        for (const x of this.elements) {
            total += x * x;
        }
        return Math.sqrt(total);
    }
    /**
     * Adds the elements of another vector to this vector and returns a new
     * vector.
     * @param v - The vector to add.
     * @returns A new vector containing the sum of the elements.
     */
    add(v) {
        const output = new Array(this.length);
        for (let i = 0; i < this.length; i++) {
            output[i] = this.elements[i] + (v.elements[i]);
        }
        return new Vector(output);
    }
    /**
     * Subtracts a vector from the current vector.
     * @param v The vector to subtract.
     * @returns A new vector representing the result of the subtraction.
     */
    subtract(v) {
        const output = new Array(this.length);
        for (let i = 0; i < this.length; i++) {
            output[i] = this.elements[i] - (v.elements[i]);
        }
        return new Vector(output);
    }
    /**
     * Scales the vector by a given factor.
     * @param n The scaling factor.
     * @returns A new Vector object representing the scaled vector.
     */
    scale(n) {
        const output = new Array(this.length);
        for (let i = 0; i < this.length; i++) {
            output[i] = this.elements[i] * n;
        }
        return new Vector(output);
    }
    /**
     * Negates the vector by scaling it by -1.
     * @returns A new Vector object representing the negated vector.
     */
    negate() {
        return this.scale(-1);
    }
    /**
     * Return the Euclidean distance between this and another Vector.
     * @param v The vector to calculate the distance to.
     * @returns The distance between the two vectors.
     */
    distance(v) {
        return this.subtract(v).magnitude();
    }
    /**
     * Normalizes the vector, making it a unit vector with the same direction but
     * a magnitude of 1. If the vector has a magnitude of 0, it returns a zero
     * vector of the same length.
     * @returns The normalized vector.
     */
    normalize() {
        const m = this.magnitude();
        if (m === 0) {
            return Vector.zero(this.length);
        }
        return this.scale(1.0 / m);
    }
    /**
     * Calculates the dot product of this vector and another vector.
     * @param v - The vector to calculate the dot product with.
     * @returns The dot product of the two vectors.
     */
    dot(v) {
        let total = 0;
        for (let i = 0; i < this.length; i++) {
            total += this.elements[i] * (v.elements[i]);
        }
        return total;
    }
    /**
     * Calculates the outer product of this vector with another vector.
     * @param v The vector to calculate the outer product with.
     * @returns A matrix representing the outer product of the two vectors.
     */
    outer(v) {
        const result = [];
        for (let i = 0; i < this.length; i++) {
            result[i] = [];
            for (let j = 0; j < v.length; j++) {
                (result[i])[j] = this.elements[i] * (v.elements[j]);
            }
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix(result);
    }
    /**
     * Calculates the cross product of this vector and the given vector.
     * @param v - The vector to calculate the cross product with.
     * @returns The resulting vector.
     */
    cross(v) {
        const output = new Array(this.length);
        for (let i = 0; i < this.length; i++) {
            output[i] =
                this.elements[(i + 1) % this.length] * (v.elements[(i + 2) % this.length]) -
                    this.elements[(i + 2) % this.length] * (v.elements[(i + 1) % this.length]);
        }
        return new Vector(output);
    }
    /**
     * Calculate the skew-symmetric matrix for this [Vector].
     * @returns The skew-symmetric matrix.
     * @throws [Error] if the vector is not of length 3.
     */
    skewSymmetric() {
        if (this.length !== 3) {
            throw new Error('Skew-symmetric matrix requires a vector of length 3.');
        }
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([
            [0, -this.elements[2], this.elements[1]],
            [this.elements[2], 0, -this.elements[0]],
            [-this.elements[1], this.elements[0], 0],
        ]);
    }
    /**
     * Rotates the vector around the X-axis by the specified angle.
     * @param theta The angle in radians.
     * @returns The rotated vector.
     */
    rotX(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const output = new Array(3);
        output[0] = this.elements[0];
        output[1] = cosT * this.elements[1] + sinT * this.elements[2];
        output[2] = -sinT * this.elements[1] + cosT * this.elements[2];
        return new Vector(output);
    }
    /**
     * Rotates the vector around the Y-axis by the specified angle.
     * @param theta The angle of rotation in radians.
     * @returns A new Vector representing the rotated vector.
     */
    rotY(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const output = new Array(3);
        output[0] = cosT * this.elements[0] + -sinT * this.elements[2];
        output[1] = this.elements[1];
        output[2] = sinT * this.elements[0] + cosT * this.elements[2];
        return new Vector(output);
    }
    /**
     * Rotates the vector around the Z-axis by the specified angle.
     * @param theta The angle of rotation in radians.
     * @returns A new Vector representing the rotated vector.
     */
    rotZ(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        const output = new Array(3);
        output[0] = cosT * this.elements[0] + sinT * this.elements[1];
        output[1] = -sinT * this.elements[0] + cosT * this.elements[1];
        output[2] = this.elements[2];
        return new Vector(output);
    }
    /**
     * Calculates the angle between this vector and another vector.
     * @param v The other vector.
     * @returns The angle between the two vectors in radians.
     */
    angle(v) {
        // better than acos for small angles
        const theta = Math.atan2(this.cross(v).magnitude(), this.dot(v));
        if (isNaN(theta)) {
            return 0.0;
        }
        return theta;
    }
    /**
     * Calculates the angle between this vector and another vector in degrees.
     * @param v The other vector.
     * @returns The angle between the two vectors in degrees.
     */
    angleDegrees(v) {
        return (this.angle(v) * (180 / Math.PI));
    }
    /**
     * Determines if there is line of sight between this vector and another vector
     * within a given radius.
     * @param v - The vector to check line of sight with.
     * @param radius - The radius within which line of sight is considered.
     * @returns True if there is line of sight, false otherwise.
     */
    sight(v, radius) {
        const r1Mag2 = this.magnitude() ** 2;
        const r2Mag2 = v.magnitude() ** 2;
        const rDot = this.dot(v);
        const tMin = (r1Mag2 - rDot) / (r1Mag2 + r2Mag2 - 2.0 * rDot);
        let los = false;
        if (tMin < 0 || tMin > 1) {
            los = true;
        }
        else {
            const c = (1.0 - tMin) * r1Mag2 + rDot * tMin;
            if (c >= radius * radius) {
                los = true;
            }
        }
        return los;
    }
    /**
     * Returns the bisect vector between this vector and the given vector. The
     * bisect vector is calculated by scaling this vector's magnitude by the
     * magnitude of the given vector, adding the result to the product of scaling
     * the given vector's magnitude by this vector's magnitude, and then
     * normalizing the resulting vector.
     * @param v - The vector to calculate the bisect with.
     * @returns The bisect vector.
     */
    bisect(v) {
        return this.scale(v.magnitude()).add(v.scale(this.magnitude())).normalize();
    }
    /**
     * Joins the current vector with another vector.
     * @param v The vector to join with.
     * @returns A new vector that contains the elements of both vectors.
     */
    join(v) {
        return new Vector(this.toList().concat(v.toList()));
    }
    /**
     * Returns a new Vector containing a portion of the elements from the
     * specified start index to the specified end index
     * @param start The start index of the slice (inclusive).
     * @param end The end index of the slice (exclusive).
     * @returns A new Vector containing the sliced elements.
     */
    slice(start, end) {
        return new Vector(this.elements.slice(start, end));
    }
    /**
     * Returns a new Matrix object representing the row vector.
     * @returns The row vector as a Matrix object.
     */
    row() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix([this.toList()]);
    }
    /**
     * Returns a new Matrix object representing the column vector of this Vector.
     * @returns The column vector as a Matrix object.
     */
    column() {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix(this.toList().map((e) => [e]));
    }
    /**
     * Converts the elements at the specified index to a Vector3D object.
     * @param index - The index of the elements to convert.
     * @returns A new Vector3D object containing the converted elements.
     */
    toVector3D(index) {
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(this.elements[index], this.elements[index + 1], this.elements[index + 2]);
    }
}


}),
"./src/engine/ootk/src/operations/Vector3D.ts": 
/*!****************************************************!*\
  !*** ./src/engine/ootk/src/operations/Vector3D.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Vector3D: () => (Vector3D)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _Matrix_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Matrix.js */ "./src/engine/ootk/src/operations/Matrix.ts");
/* ESM import */var _Vector_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Vector.js */ "./src/engine/ootk/src/operations/Vector.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



// / 3-dimensional vector.
class Vector3D {
    x;
    y;
    z;
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
        // Nothing to do here.
    }
    /**
     * Create a new Vector3D object from the first three elements of a Vector
     * object.
     * @param v The Vector object to convert.
     * @returns A new Vector3D object.
     */
    static fromVector(v) {
        return new Vector3D(v.x, v.y, v.z);
    }
    // / Origin vector.
    static origin = new Vector3D(0, 0, 0);
    // / X-axis unit vector.
    static xAxis = new Vector3D(1, 0, 0);
    // / Y-axis unit vector.
    static yAxis = new Vector3D(0, 1, 0);
    // / Z-axis unit vector.
    static zAxis = new Vector3D(0, 0, 1);
    // / Negative x-axis unit vector.
    static xAxisNeg = new Vector3D(-1, 0, 0);
    // / Negative y-axis unit vector.
    static yAxisNeg = new Vector3D(0, -1, 0);
    // / Negative z-axis unit vector.
    static zAxisNeg = new Vector3D(0, 0, -1);
    // / Convert this to a [List] of doubles.
    toList() {
        return [this.x, this.y, this.z];
    }
    // / Convert this to a [Float64List] object.
    toArray() {
        return new Float64Array([this.x, this.y, this.z]);
    }
    /**
     * Return the Vector3D element at the provided index.
     * @deprecated don't do this
     * @param index The index of the element to return.
     * @returns The element at the provided index.
     */
    getElement(index) {
        switch (index) {
            case 0:
                return this.x;
            case 1:
                return this.y;
            case 2:
                return this.z;
            default:
                throw new Error(`Index ${index} outside 3D vector bounds.`);
        }
    }
    // / Convert this to a [Vector] object.
    toVector() {
        return new _Vector_js__WEBPACK_IMPORTED_MODULE_2__.Vector(this.toList());
    }
    toString(fixed = -1) {
        if (fixed < 0) {
            return `[${this.toList().join(', ')}]`;
        }
        const output = this.toList().map((e) => e.toFixed(fixed));
        return `[${output.join(', ')}]`;
    }
    // / Return the magnitude of this vector.
    magnitude() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }
    // / Return the result of adding this to another [Vector3D].
    add(v) {
        return new Vector3D((this.x + v.x), (this.y + v.y), (this.z + v.z));
    }
    // / Return the result of subtracting this and another [Vector3D].
    subtract(v) {
        return new Vector3D((this.x - v.x), (this.y - v.y), (this.z - v.z));
    }
    // / Return a copy of this [Vector3D] scaled by [n];
    scale(n) {
        return new Vector3D(this.x * n, this.y * n, this.z * n);
    }
    // / Return a copy of this [Vector3D] with the elements negated.
    negate() {
        return this.scale(-1);
    }
    /**
     * Return the Euclidean distance between this and another Vector3D.
     * @param v The other Vector3D.
     * @returns The distance between this and the other Vector3D.
     */
    distance(v) {
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.linearDistance)(this, v);
    }
    /**
     * Convert this to a unit Vector3D.
     * @returns A unit Vector3D.
     */
    normalize() {
        const m = this.magnitude();
        if (m === 0) {
            return Vector3D.origin;
        }
        return new Vector3D(this.x / m, this.y / m, this.z / m);
    }
    // Calculate the dot product of this and another [Vector3D].
    dot(v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }
    // Calculate the outer product between this and another [Vector3D].
    outer(v) {
        return new _Matrix_js__WEBPACK_IMPORTED_MODULE_1__.Matrix([
            [this.x * v.x, this.x * v.y, this.x * v.z],
            [this.y * v.x, this.y * v.y, this.y * v.z],
            [this.z * v.x, this.z * v.y, this.z * v.z],
        ]);
    }
    // Calculate the cross product of this and another [Vector3D].
    cross(v) {
        return new Vector3D((this.y * v.z - this.z * v.y), (this.z * v.x - this.x * v.z), (this.x * v.y - this.y * v.x));
    }
    // Calculate the skew-symmetric matrix for this [Vector3D].
    skewSymmetric() {
        return new _Matrix_js__WEBPACK_IMPORTED_MODULE_1__.Matrix([
            [0, -this.z, this.y],
            [this.z, 0, -this.x],
            [-this.y, this.x, 0],
        ]);
    }
    /*
     * Create a copy of this [Vector3D] rotated in the x-axis by angle [theta]
     * _(rad)_.
     */
    rotX(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        return new Vector3D(this.x, cosT * this.y + sinT * this.z, -sinT * this.y + cosT * this.z);
    }
    /*
     * Create a copy of this [Vector3D] rotated in the y-axis by angle [theta]
     * _(rad)_.
     */
    rotY(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        return new Vector3D((cosT * this.x + -sinT * this.z), this.y, (sinT * this.x + cosT * this.z));
    }
    /*
     * Create a copy of this [Vector3D] rotated in the z-axis by angle [theta]
     * _(rad)_.
     */
    rotZ(theta) {
        const cosT = Math.cos(theta);
        const sinT = Math.sin(theta);
        return new Vector3D((cosT * this.x + sinT * this.y), (-sinT * this.x + cosT * this.y), this.z);
    }
    // Calculate the angle _(rad)_ between this and another [Vector3D].
    angle(v) {
        const theta = Math.atan2(this.cross(v).magnitude(), this.dot(v));
        return isNaN(theta) ? 0 : theta;
    }
    // Calculate the angle _(°)_ between this and another [Vector3D].
    angleDegrees(v) {
        return this.angle(v) * (180 / Math.PI);
    }
    /*
     * Return `true` if line-of-sight exists between this and another [Vector3D]
     * with a central body of the given [radius].
     */
    sight(v, radius) {
        const r1Mag2 = this.magnitude() ** 2;
        const r2Mag2 = v.magnitude() ** 2;
        const rDot = this.dot(v);
        const tMin = (r1Mag2 - rDot) / (r1Mag2 + r2Mag2 - 2 * rDot);
        let los = false;
        if (tMin < 0 || tMin > 1) {
            los = true;
        }
        else {
            const c = (1 - tMin) * r1Mag2 + rDot * tMin;
            if (c >= radius * radius) {
                los = true;
            }
        }
        return los;
    }
    // / Return the unit vector that bisects this and another [Vector3D].
    bisect(v) {
        return this.scale(v.magnitude()).add(v.scale(this.magnitude())).normalize();
    }
    // / Convert this [Vector3D] into a row [Matrix].
    row() {
        return new _Matrix_js__WEBPACK_IMPORTED_MODULE_1__.Matrix([[this.x, this.y, this.z]]);
    }
    // / Convert this [Vector3D] into a column [Matrix].
    column() {
        return new _Matrix_js__WEBPACK_IMPORTED_MODULE_1__.Matrix([[this.x], [this.y], [this.z]]);
    }
    // / Join this and another [Vector3D] into a new [Vector] object.
    join(v) {
        const output = new Float64Array(6);
        output[0] = this.x;
        output[1] = this.y;
        output[2] = this.z;
        output[3] = v.x;
        output[4] = v.y;
        output[5] = v.z;
        return new _Vector_js__WEBPACK_IMPORTED_MODULE_2__.Vector(output);
    }
}


}),
"./src/engine/ootk/src/operations/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/operations/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BoxMuller: () => (/* reexport safe */ _BoxMuller_js__WEBPACK_IMPORTED_MODULE_0__.BoxMuller),
  EulerAngles: () => (/* reexport safe */ _EulerAngles_js__WEBPACK_IMPORTED_MODULE_1__.EulerAngles),
  Quaternion: () => (/* reexport safe */ _Quaternion_js__WEBPACK_IMPORTED_MODULE_2__.Quaternion),
  RandomGaussianSource: () => (/* reexport safe */ _RandomGaussianSource_js__WEBPACK_IMPORTED_MODULE_3__.RandomGaussianSource)
});
/* ESM import */var _BoxMuller_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./BoxMuller.js */ "./src/engine/ootk/src/operations/BoxMuller.ts");
/* ESM import */var _EulerAngles_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./EulerAngles.js */ "./src/engine/ootk/src/operations/EulerAngles.ts");
/* ESM import */var _Quaternion_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Quaternion.js */ "./src/engine/ootk/src/operations/Quaternion.ts");
/* ESM import */var _RandomGaussianSource_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./RandomGaussianSource.js */ "./src/engine/ootk/src/operations/RandomGaussianSource.ts");






}),
"./src/engine/ootk/src/operations/operations.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/operations/operations.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Matrix: () => (/* reexport safe */ _Matrix_js__WEBPACK_IMPORTED_MODULE_0__.Matrix),
  Random: () => (/* reexport safe */ _Random_js__WEBPACK_IMPORTED_MODULE_1__.Random),
  Vector: () => (/* reexport safe */ _Vector_js__WEBPACK_IMPORTED_MODULE_2__.Vector),
  Vector3D: () => (/* reexport safe */ _Vector3D_js__WEBPACK_IMPORTED_MODULE_3__.Vector3D)
});
/* ESM import */var _Matrix_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Matrix.js */ "./src/engine/ootk/src/operations/Matrix.ts");
/* ESM import */var _Random_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./Random.js */ "./src/engine/ootk/src/operations/Random.ts");
/* ESM import */var _Vector_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Vector.js */ "./src/engine/ootk/src/operations/Vector.ts");
/* ESM import */var _Vector3D_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Vector3D.js */ "./src/engine/ootk/src/operations/Vector3D.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






}),
"./src/engine/ootk/src/optimize/DownhillSimplex.ts": 
/*!*********************************************************!*\
  !*** ./src/engine/ootk/src/optimize/DownhillSimplex.ts ***!
  \*********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DownhillSimplex: () => (DownhillSimplex)
});
/* ESM import */var _SimplexEntry_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./SimplexEntry.js */ "./src/engine/ootk/src/optimize/SimplexEntry.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Derivative-free Nelder-Mead simplex optimizer.
class DownhillSimplex {
    constructor() {
        // disable constructor
    }
    /**
     * Compute the centroid from a list of [SimplexEntry] objects, using cost
     * function [f].
     * @param f Cost function
     * @param xss Simplex entries
     * @returns The centroid.
     */
    static _centroid(f, xss) {
        const n = xss[0].points.length;
        const m = xss.length - 1;
        const output = new Float64Array(n);
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < n; j++) {
                output[j] += xss[i].points[j];
            }
        }
        for (let i = 0; i < n; i++) {
            output[i] /= m;
        }
        return new _SimplexEntry_js__WEBPACK_IMPORTED_MODULE_0__.SimplexEntry(f, output);
    }
    static _shrink(s, xss) {
        const x1 = xss[0];
        for (let i = 1; i < xss.length; i++) {
            const xi = xss[i];
            xss[i] = x1.modify(s, xi, x1);
        }
    }
    /**
     * Generate a new simplex from initial guess [x0], and an optional
     * simplex [step] value.
     * @param x0 Initial guess
     * @param step Simplex step
     * @returns The simplex.
     */
    static generateSimplex(x0, step = 0.01) {
        const output = [x0.slice(0)];
        for (let i = 0; i < x0.length; i++) {
            const tmp = x0.slice(0);
            tmp[i] += tmp[i] * step;
            output.push(tmp);
        }
        return output;
    }
    /**
     * Perform derivative-free Nelder-Mead simplex optimization to minimize the
     * cost function [f] for the initial simplex [xs].
     *
     * Optional arguments:
     *  - `xTolerance`: centroid delta termination criteria
     * - `fTolerance`: cost function delta termination criteria
     * - `maxIter`: maximum number of optimization iterations
     * - `adaptive`: use adaptive coefficients if possible
     * - `printIter`: print a debug statement after each iteration
     * @param f Cost function
     * @param xs Initial simplex
     * @param root0 Root0
     * @param root0.xTolerance Root0.xTolerance
     * @param root0.fTolerance Root0.fTolerance
     * @param root0.maxIter Root0.maxIter
     * @param root0.adaptive Root0.adaptive
     * @param root0.printIter Root0.printIter
     * @returns The optimal input value.
     */
    static solveSimplex(f, xs, { xTolerance = 1e-12, fTolerance = 1e-12, maxIter = 10000, adaptive = false, printIter = false, }) {
        let a;
        let g;
        let p;
        let s;
        const n = xs.length - 1;
        if (adaptive && n >= 2) {
            a = 1.0;
            g = 1.0 + 2.0 / n;
            p = 0.75 - 1.0 / (2.0 * n);
            s = 1.0 - 1.0 / n;
        }
        else {
            a = 1.0;
            g = 2.0;
            p = 0.5;
            s = 0.5;
        }
        let iter = 0;
        let action = 'init';
        const ordered = [];
        for (const x of xs) {
            ordered.push(new _SimplexEntry_js__WEBPACK_IMPORTED_MODULE_0__.SimplexEntry(f, x));
        }
        // eslint-disable-next-line no-constant-condition
        while (true) {
            ordered.sort((x, y) => x.score - y.score);
            const x0 = DownhillSimplex._centroid(f, ordered);
            // update exit criterea
            let xd = 0.0;
            let fd = 0.0;
            for (let i = 1; i < ordered.length; i++) {
                xd = Math.max(xd, x0.distance(ordered[i]));
                fd = Math.max(fd, Math.abs(x0.score - ordered[i].score));
            }
            if (printIter) {
                // eslint-disable-next-line no-console
                console.log(`${iter}: score=${x0.score} xd=${xd} fd=${fd} [${action}]`);
            }
            if (iter !== 0 && (xd < xTolerance || fd < fTolerance)) {
                return ordered[0].points;
            }
            if (iter >= maxIter) {
                return ordered[0].points;
            }
            iter++;
            // reflection
            const xr = x0.modify(a, x0, ordered[ordered.length - 1]);
            if (ordered[0].score <= xr.score && xr.score < ordered[ordered.length - 2].score) {
                ordered[ordered.length - 1] = xr;
                action = 'reflect';
                continue;
            }
            // expansion
            if (xr.score < ordered[0].score) {
                const xe = x0.modify(g, xr, x0);
                if (xe.score < xr.score) {
                    ordered[ordered.length - 1] = xe;
                }
                else {
                    ordered[ordered.length - 1] = xr;
                }
                action = 'expand';
                continue;
            }
            // contraction
            if (xr.score < ordered[ordered.length - 1].score) {
                const xc = x0.modify(p, xr, x0);
                if (xc.score < xr.score) {
                    ordered[ordered.length - 1] = xc;
                    action = 'contract';
                    continue;
                }
                else {
                    DownhillSimplex._shrink(s, ordered);
                    action = 'shrink';
                    continue;
                }
            }
            else if (xr.score >= ordered[ordered.length - 1].score) {
                const xc = x0.modify(p, ordered[ordered.length - 1], x0);
                if (xc.score < ordered[ordered.length - 1].score) {
                    ordered[ordered.length - 1] = xc;
                    action = 'contract';
                    continue;
                }
                else {
                    DownhillSimplex._shrink(s, ordered);
                    action = 'shrink';
                    continue;
                }
            }
        }
    }
}


}),
"./src/engine/ootk/src/optimize/GoldenSection.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/optimize/GoldenSection.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  GoldenSection: () => (GoldenSection)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Golden Section bounded single value optimizer.
class GoldenSection {
    static _grInv = 1.0 / (0.5 * (Math.sqrt(5) + 1));
    static check_(fc, fd, solveMax) {
        return solveMax ? fc > fd : fc < fd;
    }
    /**
     * Search for an optimal input value for function [f] that minimizes the
     * output value.
     *
     * Takes [lower] and [upper] input search bounds, and an optional
     * search [tolerance].
     * @param f Function to optimize
     * @param lower Lower bound
     * @param upper Upper bound
     * @param root0 Root0
     * @param root0.tolerance Root0.tolerance
     * @param root0.solveMax Root0.solveMax
     * @returns The optimal input value.
     */
    static search(f, lower, upper, { tolerance = 1e-5, solveMax = false, }) {
        let a = lower;
        let b = upper;
        let c = b - (b - a) * GoldenSection._grInv;
        let d = a + (b - a) * GoldenSection._grInv;
        while (Math.abs(b - a) > tolerance) {
            if (GoldenSection.check_(f(c), f(d), solveMax)) {
                b = d;
            }
            else {
                a = c;
            }
            c = b - (b - a) * GoldenSection._grInv;
            d = a + (b - a) * GoldenSection._grInv;
        }
        return 0.5 * (b + a);
    }
}


}),
"./src/engine/ootk/src/optimize/SimplexEntry.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/optimize/SimplexEntry.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  SimplexEntry: () => (SimplexEntry)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

class SimplexEntry {
    f_;
    points;
    score;
    x_;
    constructor(f_, points) {
        this.f_ = f_;
        this.points = points;
        this.x_ = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(points);
        this.score = this.f_(points);
    }
    getPoints() {
        return this.x_.toArray();
    }
    getScore() {
        return this.score;
    }
    modify(n, xa, xb) {
        return new SimplexEntry(this.f_, this.x_.add(xa.x_.subtract(xb.x_).scale(n)).toArray());
    }
    distance(se) {
        return this.x_.distance(se.x_);
    }
}


}),
"./src/engine/ootk/src/orbit_determination/BatchLeastSquaresOD.ts": 
/*!************************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/BatchLeastSquaresOD.ts ***!
  \************************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BatchLeastSquaresOD: () => (BatchLeastSquaresOD)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _observation_PropagatorPairs_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../observation/PropagatorPairs.js */ "./src/engine/ootk/src/observation/PropagatorPairs.ts");
/* ESM import */var _covariance_StateCovariance_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./../covariance/StateCovariance.js */ "./src/engine/ootk/src/covariance/StateCovariance.ts");
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _propagator_KeplerPropagator_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./../propagator/KeplerPropagator.js */ "./src/engine/ootk/src/propagator/KeplerPropagator.ts");
/* ESM import */var _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./../propagator/RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/* ESM import */var _BatchLeastSquaresResult_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./BatchLeastSquaresResult.js */ "./src/engine/ootk/src/orbit_determination/BatchLeastSquaresResult.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */







/**
 * Batch least squares orbit determination.
 */
class BatchLeastSquaresOD {
    observations_;
    apriori_;
    forceModel_;
    posStep_;
    velStep_;
    fastDerivatives_;
    /** Propagator pair cache, for generating observation Jacobians. */
    propPairs_;
    /**  Nominal state propagator. */
    propagator_;
    /**  State estimate during solve. */
    nominal_;
    /**  Solve start epoch. */
    start_;
    /**
     * Create a new [BatchLeastSquaresOD] object from a list of [Observation]
     * objects, an [apriori] state estimate, and an optional
     * spacecraft [forceModel].
     * @param observations_ List of observations.
     * @param apriori_ Apriori state estimate.
     * @param forceModel_ Spacecraft force model.
     * @param posStep_ Position step size.
     * @param velStep_ Velocity step size.
     * @param fastDerivatives_ Use fast derivatives.
     * @returns [BatchLeastSquaresOD] object.
     */
    constructor(observations_, apriori_, forceModel_, posStep_ = 1e-5, velStep_ = 1e-5, fastDerivatives_ = false) {
        this.observations_ = observations_;
        this.apriori_ = apriori_;
        this.forceModel_ = forceModel_;
        this.posStep_ = posStep_;
        this.velStep_ = velStep_;
        this.fastDerivatives_ = fastDerivatives_;
        this.observations_.sort((a, b) => a.epoch.posix - b.epoch.posix);
        this.start_ = this.observations_[0].epoch;
        this.propPairs_ = new _observation_PropagatorPairs_js__WEBPACK_IMPORTED_MODULE_1__.PropagatorPairs(this.posStep_, this.velStep_);
        this.forceModel_ ??= new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_3__.ForceModel().setGravity();
        this.propagator_ = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_5__.RungeKutta89Propagator(this.apriori_, this.forceModel_);
        this.nominal_ = this.propagator_.propagate(this.start_);
    }
    buildPropagator_(x0, simple) {
        const state = new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(this.nominal_.epoch, new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(x0[0], x0[1], x0[2]), new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector3D(x0[3], x0[4], x0[5]));
        if (simple) {
            return new _propagator_KeplerPropagator_js__WEBPACK_IMPORTED_MODULE_4__.KeplerPropagator(state.toClassicalElements());
        }
        return new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_5__.RungeKutta89Propagator(state, this.forceModel_);
    }
    static stateToX0_(state) {
        return (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.concat)(state.position.toArray(), state.velocity.toArray());
    }
    setPropagatorPairs_(x0) {
        const pl = this.buildPropagator_(x0, this.fastDerivatives_);
        for (let i = 0; i < 6; i++) {
            const step = this.propPairs_.step(i);
            const xh = x0.slice();
            xh[i] += step;
            const ph = this.buildPropagator_(xh, this.fastDerivatives_);
            this.propPairs_.set(i, ph, pl);
        }
    }
    /**
     * Attempt to solve a state estimate with the given root-mean-squared delta
     * [tolerance].
     * @param root0 Root initial guess.
     * @param root0.tolerance Root-mean-squared delta tolerance.
     * @param root0.maxIter Maximum number of iterations.
     * @param root0.printIter Print iterations.
     * @returns [BatchLeastSquaresResult] object.
     */
    solve({ tolerance = 1e-3, maxIter = 250, printIter = false, } = {}) {
        let breakFlag = false;
        const xNom = BatchLeastSquaresOD.stateToX0_(this.nominal_);
        let weightedRms = Infinity;
        const atwaMatInit = _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix.zero(6, 6);
        const atwbMatInit = _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix.zero(6, 1);
        let atwaMat = atwaMatInit;
        for (let iter = 0; iter < maxIter; iter++) {
            atwaMat = atwaMatInit;
            let atwbMat = atwbMatInit;
            this.propagator_ = this.buildPropagator_(xNom, false);
            this.setPropagatorPairs_(xNom);
            let rmsTotal = 0.0;
            let measCount = 0;
            for (const ob of this.observations_) {
                const noise = ob.noise;
                const aMat = ob.jacobian(this.propPairs_);
                const aMatTN = aMat.transpose().multiply(noise);
                const bMat = ob.residual(this.propagator_);
                atwaMat = atwaMat.add(aMatTN.multiply(aMat));
                atwbMat = atwbMat.add(aMatTN.multiply(bMat));
                rmsTotal += bMat.transpose().multiply(noise).multiply(bMat).elements[0][0];
                measCount += noise.rows;
            }
            const newWeightedRms = Math.sqrt(rmsTotal / measCount);
            if (printIter) {
                // eslint-disable-next-line no-console
                console.log(`${iter + 1}: rms=${newWeightedRms} x=${new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(xNom)}`);
            }
            if (Math.abs((weightedRms - newWeightedRms) / weightedRms) <= tolerance) {
                breakFlag = true;
            }
            weightedRms = newWeightedRms;
            const dX = atwaMat.inverse().multiply(atwbMat);
            for (let i = 0; i < 6; i++) {
                xNom[i] += dX.elements[i][0];
            }
            if (breakFlag) {
                break;
            }
        }
        const p = atwaMat.inverse();
        const covariance = new _covariance_StateCovariance_js__WEBPACK_IMPORTED_MODULE_2__.StateCovariance(p, _covariance_StateCovariance_js__WEBPACK_IMPORTED_MODULE_2__.CovarianceFrame.ECI);
        return new _BatchLeastSquaresResult_js__WEBPACK_IMPORTED_MODULE_6__.BatchLeastSquaresResult(this.buildPropagator_(xNom, false).propagate(this.start_), covariance, weightedRms);
    }
}


}),
"./src/engine/ootk/src/orbit_determination/BatchLeastSquaresResult.ts": 
/*!****************************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/BatchLeastSquaresResult.ts ***!
  \****************************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BatchLeastSquaresResult: () => (BatchLeastSquaresResult)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Batch least squares orbit determination result.
class BatchLeastSquaresResult {
    state;
    covariance;
    rms;
    /**
     * Create a new [BatchLeastSquaresResult] object, containing the solved
     * [state], [covariance], and root-mean-squared error [rms].
     * @param state The solved state.
     * @param covariance The solved covariance.
     * @param rms The root-mean-squared error.
     */
    constructor(state, covariance, rms) {
        this.state = state;
        this.covariance = covariance;
        this.rms = rms;
        // Nothing to do here.
    }
}


}),
"./src/engine/ootk/src/orbit_determination/GibbsIOD.ts": 
/*!*************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/GibbsIOD.ts ***!
  \*************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  GibbsIOD: () => (GibbsIOD)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./../propagator/RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



/**
 * Gibbs 3-position inital orbit determination.
 */
class GibbsIOD {
    mu;
    constructor(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        this.mu = mu;
        // Nothing to do here.
    }
    /** Abort solve if position plane exceeds this value. */
    static coplanarThreshold_ = (5.0 * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    /**
     * Attempt to create a state estimate from three inertial position vectors.
     *
     * Throws an error if the positions are not coplanar.
     * @param r1 Position vector 1.
     * @param r2 Position vector 2.
     * @param r3 Position vector 3.
     * @param t2 Time of position 2.
     * @param t3 Time of position 3.
     * @returns State estimate at time t2.
     */
    solve(r1, r2, r3, t2, t3) {
        const num = r1.normalize().dot(r2.normalize().cross(r3.normalize()));
        const alpha = _main_js__WEBPACK_IMPORTED_MODULE_0__.halfPi - Math.acos(num);
        if (Math.abs(alpha) > GibbsIOD.coplanarThreshold_) {
            throw new Error('Orbits are not coplanar.');
        }
        const r1m = r1.magnitude();
        const r2m = r2.magnitude();
        const r3m = r3.magnitude();
        const d = r1.cross(r2).add(r2.cross(r3).add(r3.cross(r1)));
        const n = r2.cross(r3).scale(r1m).add(r3.cross(r1).scale(r2m)).add(r1.cross(r2).scale(r3m));
        const b = d.cross(r2);
        const s = r1.scale(r2m - r3m).add(r2.scale(r3m - r1m).add(r3.scale(r1m - r2m)));
        const nm = n.magnitude();
        const dm = d.magnitude();
        const vm = Math.sqrt(this.mu / (nm * dm));
        const vlEci = b.scale((vm / r2m)).add(s.scale(vm));
        const pv = new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(t2, r2, vlEci);
        const forceModel = new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_1__.ForceModel().setGravity(this.mu);
        const orbit = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(pv, forceModel);
        const pv2 = new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(t2, r2, vlEci.negate());
        const orbit2 = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(pv2, forceModel);
        const estP3 = orbit.propagate(t3).position;
        const dist = estP3.subtract(r3).magnitude();
        const estP3b = orbit2.propagate(t3).position;
        const dist2 = estP3b.subtract(r3).magnitude();
        if (dist <= dist2) {
            orbit.reset();
            return orbit.propagate(t2);
        }
        orbit2.reset();
        return orbit2.propagate(t2);
    }
}


}),
"./src/engine/ootk/src/orbit_determination/GoodingIOD.ts": 
/*!***************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/GoodingIOD.ts ***!
  \***************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  GoodingIOD: () => (GoodingIOD)
});
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../propagator/RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/* ESM import */var _GibbsIOD_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./GibbsIOD.js */ "./src/engine/ootk/src/orbit_determination/GibbsIOD.ts");
/* ESM import */var _LambertIOD_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./LambertIOD.js */ "./src/engine/ootk/src/orbit_determination/LambertIOD.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */





/**
 * Gooding angles-only initial orbit determination.
 *
 * Used for orbit determination from three optical observations.
 */
class GoodingIOD {
    _mu;
    o1_;
    o2_;
    o3_;
    vObserverPosition1_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D.origin;
    vObserverPosition2_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D.origin;
    vObserverPosition3_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D.origin;
    r_ = 0.0;
    v_ = 0.0;
    t_ = 0.0;
    r1_ = 0.0;
    r2_ = 0.0;
    r3_ = 0.0;
    rho1_ = 0.0;
    rho2_ = 0.0;
    rho3_ = 0.0;
    d1_ = 0.0;
    d3_ = 0.0;
    facFiniteDiff_ = 0.0;
    _forceModel = new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity(1.0);
    constructor(o1, o2, o3, mu = _main_js__WEBPACK_IMPORTED_MODULE_1__.Earth.mu) {
        this._mu = mu;
        this.o1_ = o1;
        this.o2_ = o2;
        this.o3_ = o3;
    }
    _getPositionOnLoS2({ e1, r01, e3, r03, t13, t12, nRev, posigrade, }) {
        const p1 = this.vObserverPosition1_.add(e1.scale(r01));
        this.r1_ = p1.magnitude();
        const p3 = this.vObserverPosition3_.add(e3.scale(r03));
        this.r3_ = p3.magnitude();
        const p13 = p1.cross(p3);
        let th = Math.atan2(p13.magnitude(), p1.dot(p3));
        if (!posigrade) {
            th = _main_js__WEBPACK_IMPORTED_MODULE_1__.TAU - th;
        }
        const v1 = new Float64Array(2);
        const exitflag = _LambertIOD_js__WEBPACK_IMPORTED_MODULE_4__.LambertIOD.solve(this.r1_, this.r3_, th, t13, nRev, v1);
        if (exitflag) {
            const pn = p1.cross(p3);
            const pt = pn.cross(p1);
            let rt = pt.magnitude();
            if (!posigrade) {
                rt = -rt;
            }
            const vel1 = p1.scale(v1[0] / this.r1_).add(pt.scale(v1[1] / rt));
            const p2 = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(this.o1_.epoch, p1, vel1), this._forceModel).propagate(this.o1_.epoch.roll(t12)).position;
            return p2;
        }
        return null;
    }
    _modifyIterate(lineOfSight1, lineOfSight3) {
        const r13 = this.vObserverPosition3_.subtract(this.vObserverPosition1_);
        this.d1_ = r13.dot(lineOfSight1);
        this.d3_ = r13.dot(lineOfSight3);
        const d2 = lineOfSight1.dot(lineOfSight3);
        const d4 = 1.0 - d2 * d2;
        this.rho1_ = Math.max((this.d1_ - this.d3_ * d2) / d4, 0.0);
        this.rho3_ = Math.max((this.d1_ * d2 - this.d3_) / d4, 0.0);
    }
    _computeDerivatives({ x, y, lineOfSight1, lineOfSight3, pin, ein, t13, t12, nrev, direction, fd, gd, }) {
        const p = pin.normalize();
        const en = ein.normalize();
        const dx = this.facFiniteDiff_ * x;
        const dy = this.facFiniteDiff_ * y;
        const cm1 = this._getPositionOnLoS2({
            e1: lineOfSight1,
            r01: x - dx,
            e3: lineOfSight3,
            r03: y,
            t13,
            t12,
            nRev: nrev,
            posigrade: direction,
        }).subtract(this.vObserverPosition2_);
        const fm1 = p.dot(cm1);
        const gm1 = en.dot(cm1);
        const cp1 = this._getPositionOnLoS2({
            e1: lineOfSight1,
            r01: x + dx,
            e3: lineOfSight3,
            r03: y,
            t13,
            t12,
            nRev: nrev,
            posigrade: direction,
        }).subtract(this.vObserverPosition2_);
        const fp1 = p.dot(cp1);
        const gp1 = en.dot(cp1);
        const fx = (fp1 - fm1) / (2.0 * dx);
        const gx = (gp1 - gm1) / (2.0 * dx);
        const cm3 = this._getPositionOnLoS2({
            e1: lineOfSight1,
            r01: x,
            e3: lineOfSight3,
            r03: y - dy,
            t13,
            t12,
            nRev: nrev,
            posigrade: direction,
        }).subtract(this.vObserverPosition2_);
        const fm3 = p.dot(cm3);
        const gm3 = en.dot(cm3);
        const cp3 = this._getPositionOnLoS2({
            e1: lineOfSight1,
            r01: x,
            e3: lineOfSight3,
            r03: y + dy,
            t13,
            t12,
            nRev: nrev,
            posigrade: direction,
        }).subtract(this.vObserverPosition2_);
        const fp3 = p.dot(cp3);
        const gp3 = en.dot(cp3);
        const fy = (fp3 - fm3) / (2.0 * dy);
        const gy = (gp3 - gm3) / (2.0 * dy);
        fd[0] = fx;
        fd[1] = fy;
        gd[0] = gx;
        gd[1] = gy;
    }
    solve(r1Init, r3Init, nRev = 0, direction = true) {
        const lineOfSight1 = this.o1_.observation.lineOfSight();
        const lineOfSight2 = this.o2_.observation.lineOfSight();
        const lineOfSight3 = this.o3_.observation.lineOfSight();
        this.r_ = Math.max(r1Init, r3Init);
        this.v_ = Math.sqrt(this._mu / this.r_);
        this.t_ = this.r_ / this.v_;
        this.vObserverPosition1_ = this.o1_.site.position.scale(1.0 / this.r_);
        this.vObserverPosition2_ = this.o2_.site.position.scale(1.0 / this.r_);
        this.vObserverPosition3_ = this.o3_.site.position.scale(1.0 / this.r_);
        const maxiter = 100;
        this._solveRangeProblem({
            rho1init: r1Init / this.r_,
            rho3init: r3Init / this.r_,
            t13: this.o3_.epoch.difference(this.o1_.epoch) / this.t_,
            t12: this.o2_.epoch.difference(this.o1_.epoch) / this.t_,
            nrev: nRev,
            direction,
            lineOfSight1,
            lineOfSight2,
            lineOfSight3,
            maxIterations: maxiter,
        });
        const gibbs = new _GibbsIOD_js__WEBPACK_IMPORTED_MODULE_3__.GibbsIOD(this._mu);
        const p1 = this.vObserverPosition1_.add(lineOfSight1.scale(this.rho1_)).scale(this.r_);
        const p2 = this.vObserverPosition2_.add(lineOfSight2.scale(this.rho2_)).scale(this.r_);
        const p3 = this.vObserverPosition3_.add(lineOfSight3.scale(this.rho3_)).scale(this.r_);
        return gibbs.solve(p1, p2, p3, this.o2_.epoch, this.o3_.epoch);
    }
    _solveRangeProblem({ rho1init, rho3init, t13, t12, nrev, direction, lineOfSight1, lineOfSight2, lineOfSight3, maxIterations, }) {
        const arbf = 1e-6;
        const cvtol = 1e-14;
        this.rho1_ = rho1init;
        this.rho3_ = rho3init;
        let iter = 0;
        let stoppingCriterion = 10.0 * cvtol;
        while (iter < maxIterations && Math.abs(stoppingCriterion) > cvtol) {
            this.facFiniteDiff_ = arbf;
            const p2 = this._getPositionOnLoS2({
                e1: lineOfSight1,
                r01: this.rho1_,
                e3: lineOfSight3,
                r03: this.rho3_,
                t13,
                t12,
                nRev: nrev,
                posigrade: direction,
            });
            if (p2 === null) {
                this._modifyIterate(lineOfSight1, lineOfSight3);
            }
            else {
                this.r2_ = p2.magnitude();
                const c = p2.subtract(this.vObserverPosition2_);
                this.rho2_ = c.magnitude();
                const cr = lineOfSight2.dot(c);
                const u = lineOfSight2.cross(c);
                const p = u.cross(lineOfSight2).normalize();
                const ent = lineOfSight2.cross(p);
                const enr = ent.magnitude();
                if (enr === 0.0) {
                    return;
                }
                const en = ent.normalize();
                const fc = p.dot(c);
                const fd = new Float64Array(2);
                const gd = new Float64Array(2);
                this._computeDerivatives({
                    x: this.rho1_,
                    y: this.rho3_,
                    lineOfSight1,
                    lineOfSight3,
                    pin: p,
                    ein: en,
                    t13,
                    t12,
                    nrev,
                    direction,
                    fd,
                    gd,
                });
                const fr1 = fd[0];
                const fr3 = fd[1];
                const gr1 = gd[0];
                const gr3 = gd[1];
                const detj = fr1 * gr3 - fr3 * gr1;
                this.d3_ = (-gr3 * fc) / detj;
                this.d1_ = (gr1 * fc) / detj;
                this.rho1_ = this.rho1_ + this.d3_;
                this.rho3_ = this.rho3_ + this.d1_;
                const den = Math.max(cr, this.r2_);
                stoppingCriterion = fc / den;
            }
            ++iter;
        }
    }
}


}),
"./src/engine/ootk/src/orbit_determination/HerrickGibbsIOD.ts": 
/*!********************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/HerrickGibbsIOD.ts ***!
  \********************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  HerrickGibbsIOD: () => (HerrickGibbsIOD)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Herrik-Gibbs 3-position initial orbit determination.
 *
 * Possibly better than regular Gibbs IOD for closely spaced position
 * vectors (less than 5°).
 */
class HerrickGibbsIOD {
    mu;
    /**
     * Create a new [HerrickGibbsIOD] object with optional
     * gravitational parameter [mu].
     * @param mu Gravitational parameter (default: Earth.mu)
     */
    constructor(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        this.mu = mu;
        // Nothing to do here.
    }
    // / Attempt to create a state estimate from three inertial position vectors.
    solve(r1, t1, r2, t2, r3, t3) {
        const dt31 = t3.difference(t1);
        const dt32 = t3.difference(t2);
        const dt21 = t2.difference(t1);
        const r1m = r1.magnitude();
        const r2m = r2.magnitude();
        const r3m = r3.magnitude();
        const vA = r1.scale(-dt32 * (1.0 / (dt21 * dt31) + this.mu / (12.0 * r1m * r1m * r1m)));
        const vB = r2.scale((dt32 - dt21) * (1.0 / (dt21 * dt32) + this.mu / (12.0 * r2m * r2m * r2m)));
        const vC = r3.scale(dt21 * (1.0 / (dt32 * dt31) + this.mu / (12.0 * r3m * r3m * r3m)));
        const v2 = vA.add(vB).add(vC);
        return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(t2, r2, v2);
    }
}


}),
"./src/engine/ootk/src/orbit_determination/LambertIOD.ts": 
/*!***************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/LambertIOD.ts ***!
  \***************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  LambertIOD: () => (LambertIOD)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Lambert two-position and time initial orbit determination.
class LambertIOD {
    mu;
    constructor(mu = _main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.mu) {
        this.mu = mu;
        // Nothing to do here.
    }
    /**
     * Try to guess the short path argument given an [interceptor] and
     * [target] state.
     * @param interceptor Interceptor
     * @param target Target
     * @returns True if the short path should be used, false otherwise.
     */
    static useShortPath(interceptor, target) {
        const transN = interceptor.position.cross(target.position);
        const h = interceptor.position.cross(interceptor.velocity);
        return h.dot(transN) >= 0;
    }
    static timeOfFlight_(x, longway, mrev, minSma, speri, chord) {
        const a = minSma / (1.0 - x * x);
        let tof;
        if (Math.abs(x) < 1) {
            const beta = longway * 2.0 * Math.asin(Math.sqrt((speri - chord) / (2.0 * a)));
            const alpha = 2.0 * Math.acos(x);
            tof = a * Math.sqrt(a) * (alpha - Math.sin(alpha) - (beta - Math.sin(beta)) + 2.0 * Math.PI * mrev);
        }
        else {
            const alpha = 2.0 * Math.acosh(x);
            const beta = longway * 2.0 * Math.asinh(Math.sqrt((speri - chord) / (-2.0 * a)));
            tof = -a * Math.sqrt(-a) * (Math.sinh(alpha) - alpha - (Math.sinh(beta) - beta));
        }
        return tof;
    }
    /**
     * Attempt to solve output velocity [v1] _(km/s)_ given radii [r1] and
     * [r2] _(canonical)_, sweep angle [dth] _(rad)_, time of flight [tau]
     * _(canonical)_, and number of revolutions _(mRev)_.
     * @param r1 Radius 1
     * @param r2 Radius 2
     * @param dth Sweep angle
     * @param tau Time of flight
     * @param mRev Number of revolutions
     * @param v1 Output velocity
     * @returns True if successful, false otherwise.
     */
    static solve(r1, r2, dth, tau, mRev, v1) {
        const leftBranch = dth < Math.PI;
        let longway = 1;
        if (dth > Math.PI) {
            longway = -1;
        }
        const m = Math.abs(mRev);
        const rtof = Math.abs(tau);
        const theta = dth;
        const chord = Math.sqrt(r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * Math.cos(theta));
        const speri = 0.5 * (r1 + r2 + chord);
        const minSma = 0.5 * speri;
        const lambda = longway * Math.sqrt(1.0 - chord / speri);
        const logt = Math.log(rtof);
        let in1;
        let in2;
        let x1;
        let x2;
        if (m === 0) {
            in1 = -0.6523333;
            in2 = 0.6523333;
            x1 = Math.log(1.0 + in1);
            x2 = Math.log(1.0 + in2);
        }
        else {
            if (!leftBranch) {
                in1 = -0.523334;
                in2 = -0.223334;
            }
            else {
                in1 = 0.723334;
                in2 = 0.523334;
            }
            x1 = Math.tan((in1 * Math.PI) / 2);
            x2 = Math.tan((in2 * Math.PI) / 2);
        }
        const tof1 = LambertIOD.timeOfFlight_(in1, longway, m, minSma, speri, chord);
        const tof2 = LambertIOD.timeOfFlight_(in2, longway, m, minSma, speri, chord);
        let y1;
        let y2;
        if (m === 0) {
            y1 = Math.log(tof1) - logt;
            y2 = Math.log(tof2) - logt;
        }
        else {
            y1 = tof1 - rtof;
            y2 = tof2 - rtof;
        }
        let err = 1e20;
        let iterations = 0;
        const tol = 1e-13;
        const maxiter = 50;
        let xnew = 0.0;
        while (err > tol && iterations < maxiter) {
            xnew = (x1 * y2 - y1 * x2) / (y2 - y1);
            let xt;
            if (m === 0) {
                xt = Math.exp(xnew) - 1.0;
            }
            else {
                xt = (Math.atan(xnew) * 2.0) / Math.PI;
            }
            const tof = LambertIOD.timeOfFlight_(xt, longway, m, minSma, speri, chord);
            let ynew;
            if (m === 0) {
                ynew = Math.log(tof) - logt;
            }
            else {
                ynew = tof - rtof;
            }
            x1 = x2;
            x2 = xnew;
            y1 = y2;
            y2 = ynew;
            err = Math.abs(x1 - xnew);
            ++iterations;
        }
        if (err > tol) {
            return false;
        }
        let x;
        if (m === 0) {
            x = Math.exp(xnew) - 1.0;
        }
        else {
            x = (Math.atan(xnew) * 2.0) / Math.PI;
        }
        const sma = minSma / (1.0 - x * x);
        let eta;
        if (x < 1) {
            const alfa = 2.0 * Math.acos(x);
            const beta = longway * 2.0 * Math.asin(Math.sqrt((speri - chord) / (2.0 * sma)));
            const psi = (alfa - beta) / 2.0;
            const sinPsi = Math.sin(psi);
            const etaSq = (2.0 * sma * sinPsi * sinPsi) / speri;
            eta = Math.sqrt(etaSq);
        }
        else {
            const gamma = 2.0 * Math.acosh(x);
            const delta = longway * 2.0 * Math.asinh(Math.sqrt((chord - speri) / (2.0 * sma)));
            const psi = (gamma - delta) / 2.0;
            const sinhPsi = Math.sinh(psi);
            const etaSq = (-2.0 * sma * sinhPsi * sinhPsi) / speri;
            eta = Math.sqrt(etaSq);
        }
        const vr1 = (1.0 / eta) * Math.sqrt(1.0 / minSma) * ((2.0 * lambda * minSma) / r1 - (lambda + x * eta));
        const vt1 = (1.0 / eta) * Math.sqrt(1.0 / minSma) * Math.sqrt(r2 / r1) * Math.sin(dth / 2.0);
        v1[0] = vr1;
        v1[1] = vt1;
        return true;
    }
    /**
     * Estimate a state vector for inertial position [p1] _(km)_ given the
     * two epoch and positions.
     * @param p1 Position vector 1
     * @param p2 Position vector 2
     * @param t1 Epoch 1
     * @param t2 Epoch 2
     * @param root0 Optional parameters
     * @param root0.posigrade If true, use the positive root (default: true)
     * @param root0.nRev Number of revolutions (default: 0)
     * @returns A [J2000] object with the estimated state vector.
     */
    estimate(p1, p2, t1, t2, { posigrade = true, nRev = 0 } = {}) {
        const r1 = p1.magnitude();
        const r2 = p2.magnitude();
        const tof = t2.difference(t1);
        const r = Math.max(r1, r2);
        const v = Math.sqrt(this.mu / r);
        const t = r / v;
        let dth = p1.angle(p2);
        if (!posigrade) {
            dth = 2 * Math.PI - dth;
        }
        const vDep = new Float64Array(2);
        const exitFlag = LambertIOD.solve(r1 / r, r2 / r, dth, tof / t, nRev, vDep);
        if (exitFlag) {
            const pn = p1.cross(p2);
            const pt = pn.cross(p1);
            let rt = pt.magnitude();
            if (!posigrade) {
                rt = -rt;
            }
            const vel1 = p1.scale((v * vDep[0]) / r1).add(pt.scale((v * vDep[1]) / rt));
            return new _main_js__WEBPACK_IMPORTED_MODULE_0__.J2000(t1, p1, vel1);
        }
        return null;
    }
}


}),
"./src/engine/ootk/src/orbit_determination/ModifiedGoodingIOD.ts": 
/*!***********************************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/ModifiedGoodingIOD.ts ***!
  \***********************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  ModifiedGoodingIOD: () => (ModifiedGoodingIOD)
});
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../propagator/RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/* ESM import */var _optimize_DownhillSimplex_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./../optimize/DownhillSimplex.js */ "./src/engine/ootk/src/optimize/DownhillSimplex.ts");
/* ESM import */var _GoodingIOD_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./GoodingIOD.js */ "./src/engine/ootk/src/orbit_determination/GoodingIOD.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */





/**
 * Gooding angles-only initial orbit determination.
 *
 * Used for orbit determination from multiple optical observations.
 */
class ModifiedGoodingIOD {
    observations_;
    mu_;
    constructor(observations, mu = _main_js__WEBPACK_IMPORTED_MODULE_1__.Earth.mu) {
        this.observations_ = observations;
        this.mu_ = mu;
    }
    createInitial_(r0, rN, nRev, direction) {
        const iod = new _GoodingIOD_js__WEBPACK_IMPORTED_MODULE_4__.GoodingIOD(this.observations_[0], this.observations_[Math.floor(this.observations_.length / 2)], this.observations_[this.observations_.length - 1], this.mu_);
        return iod.solve(r0, rN, nRev, direction);
    }
    _createErrorFunction(aprioriEpoch) {
        const forceModel = new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity(this.mu_);
        const scoreFn = (x) => {
            const position = new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(x[0], x[1], x[2]);
            const velocity = new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(x[3], x[4], x[5]);
            const state = new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(aprioriEpoch, position, velocity);
            const propagator = new _propagator_RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_2__.RungeKutta89Propagator(state, forceModel);
            let total = 0;
            for (const oC of this.observations_) {
                const sC = propagator.propagate(oC.epoch);
                const pC = oC.site;
                const expected = oC.observation.lineOfSight();
                const actual = _main_js__WEBPACK_IMPORTED_MODULE_1__.RadecTopocentric.fromStateVector(sC, pC).lineOfSight();
                const error = expected.angle(actual);
                total += error;
            }
            return total;
        };
        return scoreFn;
    }
    solve(r0, rN, { nRev = 0, direction = true, posSearch = 10.0, velSearch = 0.1, tolerance = 1e-6, printIter = false, }) {
        if (this.observations_.length < 3) {
            throw new Error('At least 3 observations required for Gooding IOD.');
        }
        const init = this.createInitial_(r0, rN, nRev, direction);
        const guess = Float64Array.from([...init.position.toArray(), ...init.velocity.toArray()]);
        const solveFn = this._createErrorFunction(init.epoch);
        const simplex = [
            Float64Array.from(guess),
            Float64Array.from([guess[0] + posSearch, guess[1], guess[2], guess[3], guess[4], guess[5]]),
            Float64Array.from([guess[0], guess[1] + posSearch, guess[2], guess[3], guess[4], guess[5]]),
            Float64Array.from([guess[0], guess[1], guess[2] + posSearch, guess[3], guess[4], guess[5]]),
            Float64Array.from([guess[0], guess[1], guess[2], guess[3] + velSearch, guess[4], guess[5]]),
            Float64Array.from([guess[0], guess[1], guess[2], guess[3], guess[4] + velSearch, guess[5]]),
            Float64Array.from([guess[0], guess[1], guess[2], guess[3], guess[4], guess[5] + velSearch]),
        ];
        const result = _optimize_DownhillSimplex_js__WEBPACK_IMPORTED_MODULE_3__.DownhillSimplex.solveSimplex(solveFn, simplex, {
            adaptive: true,
            xTolerance: tolerance,
            fTolerance: tolerance,
            printIter,
        });
        return new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(init.epoch, new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(result[0], result[1], result[2]), new _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector3D(result[3], result[4], result[5]));
    }
}


}),
"./src/engine/ootk/src/orbit_determination/index.ts": 
/*!**********************************************************!*\
  !*** ./src/engine/ootk/src/orbit_determination/index.ts ***!
  \**********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  BatchLeastSquaresOD: () => (/* reexport safe */ _BatchLeastSquaresOD_js__WEBPACK_IMPORTED_MODULE_0__.BatchLeastSquaresOD),
  BatchLeastSquaresResult: () => (/* reexport safe */ _BatchLeastSquaresResult_js__WEBPACK_IMPORTED_MODULE_1__.BatchLeastSquaresResult),
  GibbsIOD: () => (/* reexport safe */ _GibbsIOD_js__WEBPACK_IMPORTED_MODULE_2__.GibbsIOD),
  GoodingIOD: () => (/* reexport safe */ _GoodingIOD_js__WEBPACK_IMPORTED_MODULE_3__.GoodingIOD),
  HerrickGibbsIOD: () => (/* reexport safe */ _HerrickGibbsIOD_js__WEBPACK_IMPORTED_MODULE_4__.HerrickGibbsIOD),
  LambertIOD: () => (/* reexport safe */ _LambertIOD_js__WEBPACK_IMPORTED_MODULE_5__.LambertIOD),
  ModifiedGoodingIOD: () => (/* reexport safe */ _ModifiedGoodingIOD_js__WEBPACK_IMPORTED_MODULE_6__.ModifiedGoodingIOD)
});
/* ESM import */var _BatchLeastSquaresOD_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./BatchLeastSquaresOD.js */ "./src/engine/ootk/src/orbit_determination/BatchLeastSquaresOD.ts");
/* ESM import */var _BatchLeastSquaresResult_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./BatchLeastSquaresResult.js */ "./src/engine/ootk/src/orbit_determination/BatchLeastSquaresResult.ts");
/* ESM import */var _GibbsIOD_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./GibbsIOD.js */ "./src/engine/ootk/src/orbit_determination/GibbsIOD.ts");
/* ESM import */var _GoodingIOD_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./GoodingIOD.js */ "./src/engine/ootk/src/orbit_determination/GoodingIOD.ts");
/* ESM import */var _HerrickGibbsIOD_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./HerrickGibbsIOD.js */ "./src/engine/ootk/src/orbit_determination/HerrickGibbsIOD.ts");
/* ESM import */var _LambertIOD_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./LambertIOD.js */ "./src/engine/ootk/src/orbit_determination/LambertIOD.ts");
/* ESM import */var _ModifiedGoodingIOD_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./ModifiedGoodingIOD.js */ "./src/engine/ootk/src/orbit_determination/ModifiedGoodingIOD.ts");









}),
"./src/engine/ootk/src/propagator/KeplerPropagator.ts": 
/*!************************************************************!*\
  !*** ./src/engine/ootk/src/propagator/KeplerPropagator.ts ***!
  \************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  KeplerPropagator: () => (KeplerPropagator)
});
/* ESM import */var _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../interpolator/VerletBlendInterpolator.js */ "./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Propagator.js */ "./src/engine/ootk/src/propagator/Propagator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



// / Kepler analytical two-body propagator.
class KeplerPropagator extends _Propagator_js__WEBPACK_IMPORTED_MODULE_2__.Propagator {
    initElements_;
    elements_;
    cacheState_;
    checkpoints_;
    constructor(initElements) {
        super();
        this.initElements_ = initElements;
        this.elements_ = initElements;
        this.cacheState_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000.fromClassicalElements(initElements);
        this.checkpoints_ = [];
    }
    get state() {
        return this.cacheState_;
    }
    propagate(epoch) {
        this.cacheState_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000.fromClassicalElements(this.elements_.propagate(epoch));
        return this.cacheState_;
    }
    reset() {
        this.elements_ = this.initElements_;
        this.cacheState_ = _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000.fromClassicalElements(this.elements_);
    }
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    maneuver(maneuver, interval = 60) {
        this.cacheState_ = maneuver.apply(this.propagate(maneuver.center));
        this.elements_ = this.cacheState_.toClassicalElements();
        return [this.cacheState_];
    }
    ephemerisManeuver(start, finish, maneuvers, interval = 60.0) {
        const tMvr = maneuvers.slice(0).filter((mvr) => mvr.center >= start || mvr.center <= finish);
        const ephemeris = [];
        if (tMvr[0].start > start) {
            ephemeris.push(this.propagate(start));
        }
        for (const mvr of tMvr) {
            while (this.cacheState_.epoch < mvr.center) {
                const step = Math.min(mvr.center.difference(this.cacheState_.epoch), interval);
                this.propagate(this.cacheState_.epoch.roll(step));
                if (this.cacheState_.epoch.posix !== mvr.center.posix) {
                    ephemeris.push(this.cacheState_);
                }
            }
            ephemeris.push(...this.maneuver(mvr, interval));
        }
        while (this.cacheState_.epoch < finish) {
            const step = Math.min(finish.difference(this.cacheState_.epoch), interval);
            this.propagate(this.cacheState_.epoch.roll(step));
            ephemeris.push(this.cacheState_);
        }
        return new _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_0__.VerletBlendInterpolator(ephemeris);
    }
    checkpoint() {
        this.checkpoints_.push(this.cacheState_);
        return this.checkpoints_.length - 1;
    }
    clearCheckpoints() {
        this.checkpoints_ = [];
    }
    restore(index) {
        this.cacheState_ = this.checkpoints_[index];
    }
}


}),
"./src/engine/ootk/src/propagator/Propagator.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/propagator/Propagator.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Propagator: () => (Propagator)
});
/* ESM import */var _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../interpolator/VerletBlendInterpolator.js */ "./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _optimize_GoldenSection_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./../optimize/GoldenSection.js */ "./src/engine/ootk/src/optimize/GoldenSection.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



// Propagator base class.
class Propagator {
    /*
     * Generate a [VerletBlendInterpolator] containing ephemeris over the
     * [start] and [stop] propagation period, with an optional
     * ephemeris [interval].
     */
    ephemeris(start, stop, interval = 60.0) {
        const output = [this.propagate(start)];
        let tempEpoch = start;
        while (tempEpoch <= stop) {
            tempEpoch = tempEpoch.roll(interval);
            output.push(this.propagate(tempEpoch));
        }
        return new _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_0__.VerletBlendInterpolator(output);
    }
    /*
     * Generate a list of [J2000] states integrating over a maneuver.
     *
     * If the maneuver is impulsive, the list will only contain a single state.
     */
    maneuver(maneuver, interval = 60.0) {
        const output = [];
        const tempEpoch = maneuver.start;
        const stop = maneuver.stop;
        output.push(this.propagate(tempEpoch));
        let current = tempEpoch;
        while (current <= stop) {
            current = current.roll(interval);
            output.push(this.propagate(current));
        }
        return output;
    }
    /*
     * Generate a [VerletBlendInterpolator] containing maneuver ephemeris over
     * the [start] and [finish] interval, with an optional ephemeris [interval].
     */
    ephemerisManeuver(start, finish, maneuvers, interval = 60.0) {
        const output = [];
        const tempEpoch = start;
        output.push(this.propagate(tempEpoch));
        const stop = finish;
        let current = tempEpoch;
        while (current <= stop) {
            current = current.roll(interval);
            output.push(this.propagate(current));
        }
        return new _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_0__.VerletBlendInterpolator(output);
    }
    // Return the epoch of the ascending node after the [start] epoch.
    ascendingNodeEpoch(start) {
        const period = this.state.period / 60;
        const step = period / 8;
        let current = start;
        const stop = current.roll(period);
        this.propagate(current);
        let previous = this.state.position.z;
        while (current <= stop) {
            current = current.roll(step);
            this.propagate(current);
            if (Math.sign(this.state.position.z) === Math.sign(-previous) && this.state.velocity.z > 0) {
                break;
            }
            previous = this.state.position.z;
        }
        const t = _optimize_GoldenSection_js__WEBPACK_IMPORTED_MODULE_2__.GoldenSection.search((x) => {
            this.propagate(new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(x));
            return Math.abs(this.state.position.z);
        }, current.posix - step, current.posix, { tolerance: 1e-3 });
        return new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(t);
    }
    // Return the epoch of the descending node after the [start] epoch.
    descendingNodeEpoch(start) {
        const period = this.state.period / 60;
        const step = period / 8;
        let current = start;
        const stop = current.roll(period);
        this.propagate(current);
        let previous = this.state.position.z;
        while (current <= stop) {
            current = current.roll(step);
            this.propagate(current);
            if (Math.sign(this.state.position.z) === Math.sign(-previous) && this.state.velocity.z < 0) {
                break;
            }
            previous = this.state.position.z;
        }
        const t = _optimize_GoldenSection_js__WEBPACK_IMPORTED_MODULE_2__.GoldenSection.search((x) => {
            this.propagate(new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(x));
            return Math.abs(this.state.position.z);
        }, current.posix - step, current.posix, { tolerance: 1e-3 });
        return new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(t);
    }
    // Return the epoch of apogee after the [start] epoch.
    apogeeEpoch(start) {
        const slice = 8;
        const period = this.state.period / 60;
        const step = period / slice;
        let current = start;
        this.propagate(current);
        let tCache = current;
        let rCache = this.state.position.magnitude();
        for (let i = 0; i < slice; i++) {
            current = current.roll(step);
            const t = new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(_optimize_GoldenSection_js__WEBPACK_IMPORTED_MODULE_2__.GoldenSection.search((x) => {
                this.propagate(new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(x));
                return this.state.position.magnitude();
            }, current.posix - step, current.posix, { tolerance: 1e-3, solveMax: true }));
            this.propagate(t);
            const r = this.state.position.magnitude();
            if (r > rCache) {
                tCache = t;
                rCache = r;
            }
        }
        return tCache;
    }
    // Return the epoch of perigee after the [start] epoch.
    perigeeEpoch(start) {
        const slice = 8;
        const period = this.state.period / 60;
        const step = period / slice;
        let current = start;
        this.propagate(current);
        let tCache = current;
        let rCache = this.state.position.magnitude();
        for (let i = 0; i < slice; i++) {
            current = current.roll(step);
            const t = new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(_optimize_GoldenSection_js__WEBPACK_IMPORTED_MODULE_2__.GoldenSection.search((x) => {
                this.propagate(new _main_js__WEBPACK_IMPORTED_MODULE_1__.EpochUTC(x));
                return this.state.position.magnitude();
            }, current.posix - step, current.posix, { tolerance: 1e-3, solveMax: false }));
            this.propagate(t);
            const r = this.state.position.magnitude();
            if (r < rCache) {
                tCache = t;
                rCache = r;
            }
        }
        return tCache;
    }
}


}),
"./src/engine/ootk/src/propagator/RkCheckpoint.ts": 
/*!********************************************************!*\
  !*** ./src/engine/ootk/src/propagator/RkCheckpoint.ts ***!
  \********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RkCheckpoint: () => (RkCheckpoint)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Runge-Kutta adaptive state checkpoint.
class RkCheckpoint {
    cacheState;
    stepSize;
    constructor(cacheState, stepSize) {
        this.cacheState = cacheState;
        this.stepSize = stepSize;
        // Nothing to do here.
    }
}


}),
"./src/engine/ootk/src/propagator/RkResult.ts": 
/*!****************************************************!*\
  !*** ./src/engine/ootk/src/propagator/RkResult.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RkResult: () => (RkResult)
});
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// / Result of adaptive numerical integration.
class RkResult {
    state;
    error;
    newStep;
    // / Create a new [RkResult] object.
    constructor(state, error, newStep) {
        this.state = state;
        this.error = error;
        this.newStep = newStep;
        // Nothing to do here.
    }
}


}),
"./src/engine/ootk/src/propagator/RungeKutta4Propagator.ts": 
/*!*****************************************************************!*\
  !*** ./src/engine/ootk/src/propagator/RungeKutta4Propagator.ts ***!
  \*****************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RungeKutta4Propagator: () => (RungeKutta4Propagator)
});
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../interpolator/VerletBlendInterpolator.js */ "./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _Propagator_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Propagator.js */ "./src/engine/ootk/src/propagator/Propagator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */




// / Runge-Kutta 4 fixed numerical propagator.
class RungeKutta4Propagator extends _Propagator_js__WEBPACK_IMPORTED_MODULE_3__.Propagator {
    initState_;
    forceModel_;
    stepSize_;
    cacheState_;
    checkpoints_;
    /**
     * Create a new [RungeKutta4Propagator] object from an initial state vector and
     * along with an optional [ForceModel] and [stepSize] in seconds.
     * @param initState_ Initial state vector.
     * @param forceModel_ Numerical integration force model.
     * @param stepSize_ Integration step size _(seconds)_.
     * @param cacheState_ Cached state vector.
     * @param checkpoints_ Cached state vector checkpoints.
     */
    constructor(initState_, forceModel_ = new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity(), stepSize_ = 15.0, cacheState_ = initState_, checkpoints_ = []) {
        super();
        this.initState_ = initState_;
        this.forceModel_ = forceModel_;
        this.stepSize_ = stepSize_;
        this.cacheState_ = cacheState_;
        this.checkpoints_ = checkpoints_;
        this.stepSize_ = Math.abs(stepSize_);
    }
    // / Set the integrator step size to the provided number of [seconds].
    setStepSize(seconds) {
        this.stepSize_ = Math.abs(seconds);
    }
    // / Set numerical integration force model.
    setForceModel(forceModel) {
        this.forceModel_ = forceModel;
    }
    ephemerisManeuver(start, finish, maneuvers, interval = 60.0) {
        const tMvr = maneuvers.slice(0).filter((mvr) => mvr.start >= start || mvr.stop <= finish);
        const ephemeris = [];
        if (tMvr[0].start > start) {
            ephemeris.push(this.propagate(start));
        }
        for (const mvr of tMvr) {
            while (this.cacheState_.epoch < mvr.start) {
                const step = Math.min(mvr.start.difference(this.cacheState_.epoch), interval);
                this.propagate(this.cacheState_.epoch.roll(step));
                if (this.cacheState_.epoch.posix !== mvr.start.posix) {
                    ephemeris.push(this.cacheState_);
                }
            }
            ephemeris.push(...this.maneuver(mvr, interval));
        }
        while (this.cacheState_.epoch.posix < finish.posix) {
            const step = Math.min(finish.difference(this.cacheState_.epoch), interval);
            this.propagate(this.cacheState_.epoch.roll(step));
            ephemeris.push(this.cacheState_);
        }
        return new _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_1__.VerletBlendInterpolator(ephemeris);
    }
    maneuver(maneuver, interval = 60.0) {
        if (maneuver.isImpulsive) {
            this.cacheState_ = maneuver.apply(this.propagate(maneuver.center));
            return [this.cacheState_];
        }
        let tState = this.propagate(maneuver.start);
        this.forceModel_.loadManeuver(maneuver);
        const ephemeris = [tState];
        while (tState.epoch < maneuver.stop) {
            const step = Math.min(maneuver.stop.difference(tState.epoch), interval);
            tState = this.propagate(tState.epoch.roll(step));
            ephemeris.push(tState);
        }
        this.forceModel_.clearManeuver();
        return ephemeris;
    }
    _kFn(state, hArg, kArg) {
        const epoch = state.epoch.roll(hArg);
        const posvel = state.position.join(state.velocity);
        const result = posvel.add(kArg);
        const sample = new _main_js__WEBPACK_IMPORTED_MODULE_2__.J2000(epoch, result.toVector3D(0), result.toVector3D(3));
        return this.forceModel_.derivative(sample);
    }
    _integrate(state, step) {
        const k1 = this._kFn(state, 0, _main_js__WEBPACK_IMPORTED_MODULE_2__.Vector.zero(6)).scale(step);
        const k2 = this._kFn(state, 0.5 * step, k1.scale(0.5)).scale(step);
        const k3 = this._kFn(state, 0.5 * step, k2.scale(0.5)).scale(step);
        const k4 = this._kFn(state, step, k3).scale(step);
        const v1 = k1;
        const v2 = v1.add(k2.scale(2));
        const v3 = v2.add(k3.scale(2));
        const v4 = v3.add(k4);
        const tNext = state.epoch.roll(step);
        const posvel = state.position.join(state.velocity);
        const result = posvel.add(v4.scale(1 / 6));
        return new _main_js__WEBPACK_IMPORTED_MODULE_2__.J2000(tNext, result.toVector3D(0), result.toVector3D(3));
    }
    propagate(epoch) {
        let delta = epoch.difference(this.cacheState_.epoch);
        while (delta !== 0) {
            const direction = delta >= 0 ? 1 : -1;
            const dt = Math.min(Math.abs(delta), this.stepSize_) * direction;
            this.cacheState_ = this._integrate(this.cacheState_, dt);
            delta = epoch.difference(this.cacheState_.epoch);
        }
        return this.cacheState_;
    }
    reset() {
        this.cacheState_ = this.initState_;
    }
    get state() {
        return this.cacheState_;
    }
    checkpoint() {
        this.checkpoints_.push(this.cacheState_);
        return this.checkpoints_.length - 1;
    }
    clearCheckpoints() {
        this.checkpoints_ = [];
    }
    restore(index) {
        this.cacheState_ = this.checkpoints_[index];
    }
}


}),
"./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts": 
/*!******************************************************************!*\
  !*** ./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts ***!
  \******************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RungeKutta89Propagator: () => (RungeKutta89Propagator)
});
/* ESM import */var _RungeKuttaAdaptive_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./RungeKuttaAdaptive.js */ "./src/engine/ootk/src/propagator/RungeKuttaAdaptive.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/* eslint-disable @typescript-eslint/no-loss-of-precision */
/* eslint-disable class-methods-use-this */

// / Runge-Kutta 8(9) adaptive numerical propagator.
class RungeKutta89Propagator extends _RungeKuttaAdaptive_js__WEBPACK_IMPORTED_MODULE_0__.RungeKuttaAdaptive {
    a_ = Float64Array.from([
        0.0, 0.3462e-1, 0.9702435063878044594828361677100617517633e-1, 0.1455365259581706689224254251565092627645, 0.561,
        0.2290079115904850126662751771814700052182, 0.5449920884095149873337248228185299947818, 0.645, 0.48375, 0.6757e-1,
        0.25, 0.6590650618730998549405331618649220295334, 0.8206, 0.9012, 1.0, 1.0,
    ]);
    b_ = [
        new Float64Array(0),
        Float64Array.from([0.3462e-1]),
        Float64Array.from([-0.389335438857287327017042687229284478532e-1, 0.1359578945245091786499878854939346230295]),
        Float64Array.from([0.3638413148954266723060635628912731569111e-1, 0.0, 0.1091523944686280016918190688673819470733]),
        Float64Array.from([
            2.025763914393969636805657604282571047511, 0.0, -7.638023836496292020387602153091964592952,
            6.173259922102322383581944548809393545442,
        ]),
        Float64Array.from([
            0.5112275589406060872792270881648288397197e-1, 0.0, 0.0, 0.1770823794555021537929910813839068684087,
            0.80277624092225014536138698108025283759e-3,
        ]),
        Float64Array.from([
            0.1316006357975216279279871693164256985334, 0.0, 0.0, -0.2957276252669636417685183174672273730699,
            0.878137803564295237421124704053886667082e-1, 0.62130529752252747743214350056394300261,
        ]),
        Float64Array.from([
            0.7166666666666666666666666666666666666667e-1, 0.0, 0.0, 0.0, 0.0, 0.3305533578915319409260346730051472207728,
            0.2427799754418013924072986603281861125606,
        ]),
        Float64Array.from([
            0.71806640625e-1, 0.0, 0.0, 0.0, 0.0, 0.3294380283228177160744825466257672816401,
            0.1165190029271822839255174533742327183599, -0.34013671875e-1,
        ]),
        Float64Array.from([
            0.4836757646340646986611287718844085773549e-1, 0.0, 0.0, 0.0, 0.0, 0.3928989925676163974333190042057047002852e-1,
            0.1054740945890344608263649267140088017604, -0.2143865284648312665982642293830533996214e-1,
            -0.1041229174627194437759832813847147895623,
        ]),
        Float64Array.from([
            -0.2664561487201478635337289243849737340534e-1, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333333333333333333333333333e-1,
            -0.1631072244872467239162704487554706387141, 0.3396081684127761199487954930015522928244e-1,
            0.1572319413814626097110769806810024118077, 0.215226747803187955230353477879477037696,
        ]),
        Float64Array.from([
            0.3689009248708622334786359863227633989718e-1, 0.0, 0.0, 0.0, 0.0, -0.1465181576725542928653609891758501156785,
            0.2242577768172024345345469822625833796001, 0.2294405717066072637090897902753790803034e-1,
            -0.35850052905728761357394424889330334334e-2, 0.8669223316444385506869203619044453906053e-1,
            0.4383840651968337846196219974168630120572,
        ]),
        Float64Array.from([
            -0.4866012215113340846662212357570395295088, 0.0, 0.0, 0.0, 0.0, -6.304602650282852990657772792012007122988,
            -0.281245618289472564778284183790118418111, -2.679019236219849057687906597489223155566,
            0.518815663924157511565311164615012522024, 1.365353187603341710683633635235238678626,
            5.885091088503946585721274891680604830712, 2.802808786272062889819965117517532194812,
        ]),
        Float64Array.from([
            0.4185367457753471441471025246471931649633, 0.0, 0.0, 0.0, 0.0, 6.724547581906459363100870806514855026676,
            -0.425444280164611790606983409697113064616, 3.343279153001265577811816947557982637749,
            0.617081663117537759528421117507709784737, -0.929966123939932833937749523988800852013,
            -6.099948804751010722472962837945508844846, -3.002206187889399044804158084895173690015,
            0.2553202529443445472336424602988558373637,
        ]),
        Float64Array.from([
            -0.779374086122884664644623040843840506343, 0.0, 0.0, 0.0, 0.0, -13.93734253810777678786523664804936051203,
            1.252048853379357320949735183924200895136, -14.69150040801686878191527989293072091588,
            -0.494705058533141685655191992136962873577, 2.242974909146236657906984549543692874755,
            13.36789380382864375813864978592679139881, 14.39665048665068644512236935340272139005,
            -0.7975813331776800379127866056663258667437, 0.4409353709534277758753793068298041158235,
        ]),
        Float64Array.from([
            2.058051337466886442151242368989994043993, 0.0, 0.0, 0.0, 0.0, 22.35793772796803295519317565842520212899,
            0.90949810997556332745009198137971890783, 35.89110098240264104710550686568482456493,
            -3.442515027624453437985000403608480262211, -4.865481358036368826566013387928704014496,
            -18.90980381354342625688427480879773032857, -34.26354448030451782929251177395134170515,
            1.264756521695642578827783499806516664686, 0.0, 0.0,
        ]),
    ];
    ch_ = Float64Array.from([
        0.1996996514886773085518508418098868756464e-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        2.191499304949330054530747099310837524864, 0.8857071848208438030833722031786358862953e-1,
        0.1140560234865965622484956605091432032674, 0.2533163805345107065564577734569651977347,
        -2.056564386240941011158999594595981300493, 0.340809679901311993516009489422454381283, 0.0, 0.0,
        0.4834231373823958314376726739772871714902e-1,
    ]);
    c_ = Float64Array.from([
        0.1461197685842315252051541915018784713459e-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        -0.391521186233133908941022826728824203081, 0.2310932500289506415909675644868993669908,
        0.1274766769992852382560589467488989175618, 0.2246434176204157731566981937082069688984,
        0.5684352689748512932705226972873692126743, 0.5825871557215827200814768021863420902155e-1,
        0.1364317403482215641609022744494239843327, 0.3057013983082797397721005067920369646664e-1, 0.0,
    ]);
    get a() {
        return this.a_;
    }
    get b() {
        return this.b_;
    }
    get ch() {
        return this.ch_;
    }
    get c() {
        return this.c_;
    }
    get order() {
        return 8;
    }
}


}),
"./src/engine/ootk/src/propagator/RungeKuttaAdaptive.ts": 
/*!**************************************************************!*\
  !*** ./src/engine/ootk/src/propagator/RungeKuttaAdaptive.ts ***!
  \**************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RungeKuttaAdaptive: () => (RungeKuttaAdaptive)
});
/* ESM import */var _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../force/ForceModel.js */ "./src/engine/ootk/src/force/ForceModel.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./../interpolator/VerletBlendInterpolator.js */ "./src/engine/ootk/src/interpolator/VerletBlendInterpolator.ts");
/* ESM import */var _Propagator_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Propagator.js */ "./src/engine/ootk/src/propagator/Propagator.ts");
/* ESM import */var _RkCheckpoint_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./RkCheckpoint.js */ "./src/engine/ootk/src/propagator/RkCheckpoint.ts");
/* ESM import */var _RkResult_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./RkResult.js */ "./src/engine/ootk/src/propagator/RkResult.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */






// / Adaptive Runge-Kutta propagator base class.
class RungeKuttaAdaptive extends _Propagator_js__WEBPACK_IMPORTED_MODULE_3__.Propagator {
    initState_;
    forceModel_;
    tolerance_;
    /**
     * Create a new [RungeKuttaAdaptive] object from an initial state vector
     * along with an optional [ForceModel] and [tolerance].
     * @param initState_ Initial state vector.
     * @param forceModel_ Numerical integration force model.
     * @param tolerance_ Minimum allowable local error tolerance.
     */
    constructor(initState_, forceModel_ = new _force_ForceModel_js__WEBPACK_IMPORTED_MODULE_0__.ForceModel().setGravity(), tolerance_ = 1e-9) {
        super();
        this.initState_ = initState_;
        this.forceModel_ = forceModel_;
        this.tolerance_ = tolerance_;
        this._cacheState = this.initState_;
        this.tolerance_ = Math.max(RungeKuttaAdaptive._minTolerance, Math.abs(tolerance_));
    }
    // / Initial state vector.
    _cacheState;
    _checkpoints = [];
    // / Integration step size _(seconds)_.
    _stepSize = 60.0;
    // / Minimum allowable local error tolerance.
    static _minTolerance = 1e-15;
    get state() {
        return this._cacheState;
    }
    reset() {
        this._cacheState = this.initState_;
        this._stepSize = 60.0;
    }
    // / Set numerical integration force model.
    setForceModel(forceModel) {
        this.forceModel_ = forceModel;
    }
    kfn_(epoch, rv, hArg, kArg, step) {
        const t = epoch.roll(hArg * step);
        const rvNew = rv.add(kArg);
        const sample = new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(t, rvNew.toVector3D(0), rvNew.toVector3D(3));
        return this.forceModel_.derivative(sample).scale(step);
    }
    integrate_(state, step) {
        const k = new Array(this.a.length).fill(_main_js__WEBPACK_IMPORTED_MODULE_1__.Vector.origin3);
        const y = state.position.join(state.velocity);
        for (let i = 0; i < this.a.length; i++) {
            let kArg = _main_js__WEBPACK_IMPORTED_MODULE_1__.Vector.origin6;
            if (i !== 0) {
                for (let j = 0; j < i; j++) {
                    kArg = kArg.add(k[j].scale(this.b[i][j]));
                }
            }
            k[i] = this.kfn_(state.epoch, y, this.a[i], kArg, step);
        }
        let y1 = y;
        let y2 = y;
        for (let i = 0; i < k.length; i++) {
            y1 = y1.add(k[i].scale(this.ch[i]));
            y2 = y2.add(k[i].scale(this.c[i]));
        }
        const teVal = y1.distance(y2);
        let hNew = 0.9 * step * (this.tolerance_ / teVal) ** (1.0 / this.order);
        const hOld = Math.abs(step);
        hNew = Math.max(0.2 * hOld, Math.min(5.0 * hOld, hNew));
        hNew = Math.max(1e-5, Math.min(1000.0, hNew));
        return new _RkResult_js__WEBPACK_IMPORTED_MODULE_5__.RkResult(new _main_js__WEBPACK_IMPORTED_MODULE_1__.J2000(state.epoch.roll(step), y1.toVector3D(0), y1.toVector3D(3)), teVal, hNew);
    }
    propagate(epoch) {
        let delta = epoch.difference(this._cacheState.epoch);
        while (delta !== 0) {
            const direction = delta >= 0 ? 1 : -1;
            const dt = Math.min(Math.abs(delta), this._stepSize) * direction;
            const result = this.integrate_(this._cacheState, dt);
            this._stepSize = result.newStep;
            if (result.error > this.tolerance_) {
                continue;
            }
            this._cacheState = result.state;
            delta = epoch.difference(this._cacheState.epoch);
        }
        return this._cacheState;
    }
    maneuver(maneuver, interval = 60.0) {
        if (maneuver.isImpulsive) {
            this._cacheState = maneuver.apply(this.propagate(maneuver.center));
            return [this._cacheState];
        }
        let tState = this.propagate(maneuver.start);
        this.forceModel_.loadManeuver(maneuver);
        const ephemeris = [tState];
        while (tState.epoch < maneuver.stop) {
            const step = Math.min(maneuver.stop.difference(tState.epoch), interval);
            tState = this.propagate(tState.epoch.roll(step));
            ephemeris.push(tState);
        }
        this.forceModel_.clearManeuver();
        return ephemeris;
    }
    ephemerisManeuver(start, finish, maneuvers, interval = 60.0) {
        const tMvr = maneuvers.filter((mvr) => mvr.start >= start || mvr.stop <= finish);
        const ephemeris = [];
        if (tMvr[0].start > start) {
            ephemeris.push(this.propagate(start));
        }
        for (const mvr of tMvr) {
            while (this._cacheState.epoch < mvr.start) {
                const step = Math.min(mvr.start.difference(this._cacheState.epoch), interval);
                this.propagate(this._cacheState.epoch.roll(step));
                if (this._cacheState.epoch.posix !== mvr.start.posix) {
                    ephemeris.push(this._cacheState);
                }
            }
            ephemeris.push(...this.maneuver(mvr, interval));
        }
        while (this._cacheState.epoch.posix < finish.posix) {
            const step = Math.min(finish.difference(this._cacheState.epoch), interval);
            this.propagate(this._cacheState.epoch.roll(step));
            ephemeris.push(this._cacheState);
        }
        return new _interpolator_VerletBlendInterpolator_js__WEBPACK_IMPORTED_MODULE_2__.VerletBlendInterpolator(ephemeris);
    }
    checkpoint() {
        this._checkpoints.push(new _RkCheckpoint_js__WEBPACK_IMPORTED_MODULE_4__.RkCheckpoint(this._cacheState, this._stepSize));
        return this._checkpoints.length - 1;
    }
    clearCheckpoints() {
        this._checkpoints.length = 0;
    }
    restore(index) {
        const checkpoint = this._checkpoints[index];
        this._cacheState = checkpoint.cacheState;
        this._stepSize = checkpoint.stepSize;
    }
}


}),
"./src/engine/ootk/src/propagator/Sgp4Propagator.ts": 
/*!**********************************************************!*\
  !*** ./src/engine/ootk/src/propagator/Sgp4Propagator.ts ***!
  \**********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sgp4Propagator: () => (Sgp4Propagator)
});
/* ESM import */var _Propagator_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Propagator.js */ "./src/engine/ootk/src/propagator/Propagator.ts");
/**
 * @author Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Sgp4Propagator is a propagator that uses the SGP4 model to propagate the state of an object.
 * This class is useful for propagating multiple states with the same TLE, since it caches the
 * state of the TLE at different epochs.
 */
class Sgp4Propagator extends _Propagator_js__WEBPACK_IMPORTED_MODULE_0__.Propagator {
    tle_;
    constructor(tle_) {
        super();
        this.tle_ = tle_;
        this.cacheState_ = tle_.state.toJ2000();
    }
    cacheState_;
    checkpoints_ = [];
    /**
     * Gets the state of the propagator in the J2000 coordinate system.
     * @returns The J2000 state of the propagator.
     */
    get state() {
        return this.cacheState_;
    }
    /**
     * Calculates the ephemeris maneuver using the SGP4 propagator.
     * @param start The start epoch in UTC.
     * @param finish The finish epoch in UTC.
     * @param maneuvers The array of thrust maneuvers.
     * @param interval The time interval in seconds.
     */
    ephemerisManeuver(start, finish, maneuvers, interval = 60.0) {
        throw new Error('Maneuvers cannot be modelled with SGP4.');
    }
    /**
     * Performs a maneuver with the given thrust.
     * @param maneuver - The thrust maneuver to perform.
     * @param interval - The time interval for the maneuver (default: 60.0 seconds).
     * @throws Error if maneuvers cannot be modeled with SGP4.
     */
    maneuver(maneuver, interval = 60.0) {
        throw new Error('Maneuvers cannot be modelled with SGP4.');
    }
    /**
     * Propagates the state of the Sgp4Propagator to a specified epoch in J2000 coordinates.
     * @param epoch - The epoch in UTC format.
     * @returns The propagated state in J2000 coordinates.
     */
    propagate(epoch) {
        this.cacheState_ = this.tle_.propagate(epoch).toJ2000();
        return this.cacheState_;
    }
    /**
     * Resets the state of the Sgp4Propagator by updating the cache state
     * to the current J2000 state of the TLE.
     */
    reset() {
        this.cacheState_ = this.tle_.state.toJ2000();
    }
    /**
     * Saves the current state of the propagator and returns the index of the checkpoint.
     * @returns The index of the checkpoint.
     */
    checkpoint() {
        this.checkpoints_.push(this.cacheState_);
        return this.checkpoints_.length - 1;
    }
    /**
     * Clears all the checkpoints in the propagator.
     */
    clearCheckpoints() {
        this.checkpoints_ = [];
    }
    /**
     * Restores the state of the propagator to a previously saved checkpoint.
     * @param index - The index of the checkpoint to restore.
     */
    restore(index) {
        this.cacheState_ = this.checkpoints_[index];
    }
}


}),
"./src/engine/ootk/src/propagator/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/propagator/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  RungeKutta4Propagator: () => (/* reexport safe */ _RungeKutta4Propagator_js__WEBPACK_IMPORTED_MODULE_0__.RungeKutta4Propagator),
  RungeKutta89Propagator: () => (/* reexport safe */ _RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_1__.RungeKutta89Propagator),
  Sgp4Propagator: () => (/* reexport safe */ _Sgp4Propagator_js__WEBPACK_IMPORTED_MODULE_2__.Sgp4Propagator)
});
/* ESM import */var _RungeKutta4Propagator_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./RungeKutta4Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta4Propagator.ts");
/* ESM import */var _RungeKutta89Propagator_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./RungeKutta89Propagator.js */ "./src/engine/ootk/src/propagator/RungeKutta89Propagator.ts");
/* ESM import */var _Sgp4Propagator_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./Sgp4Propagator.js */ "./src/engine/ootk/src/propagator/Sgp4Propagator.ts");





}),
"./src/engine/ootk/src/sgp4/index.ts": 
/*!*******************************************!*\
  !*** ./src/engine/ootk/src/sgp4/index.ts ***!
  \*******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sgp4: () => (/* reexport safe */ _sgp4_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4)
});
/* ESM import */var _sgp4_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./sgp4.js */ "./src/engine/ootk/src/sgp4/sgp4.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */



}),
"./src/engine/ootk/src/sgp4/sgp4-error.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/sgp4/sgp4-error.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sgp4Error: () => (Sgp4Error),
  Sgp4ErrorCode: () => (Sgp4ErrorCode)
});
/**
 * Improved error handling for SGP4
 */
// Define specific error types
var Sgp4ErrorCode;
(function (Sgp4ErrorCode) {
    Sgp4ErrorCode[Sgp4ErrorCode["NO_ERROR"] = 0] = "NO_ERROR";
    Sgp4ErrorCode[Sgp4ErrorCode["MEAN_ELEMENTS_INVALID"] = 1] = "MEAN_ELEMENTS_INVALID";
    Sgp4ErrorCode[Sgp4ErrorCode["MEAN_MOTION_NEGATIVE"] = 2] = "MEAN_MOTION_NEGATIVE";
    Sgp4ErrorCode[Sgp4ErrorCode["PERT_ELEMENTS_INVALID"] = 3] = "PERT_ELEMENTS_INVALID";
    Sgp4ErrorCode[Sgp4ErrorCode["SEMI_LATUS_RECTUM_NEGATIVE"] = 4] = "SEMI_LATUS_RECTUM_NEGATIVE";
    Sgp4ErrorCode[Sgp4ErrorCode["EPOCH_ELEMENTS_SUBORBITAL"] = 5] = "EPOCH_ELEMENTS_SUBORBITAL";
    Sgp4ErrorCode[Sgp4ErrorCode["SATELLITE_DECAYED"] = 6] = "SATELLITE_DECAYED"; // satellite has decayed
})(Sgp4ErrorCode || (Sgp4ErrorCode = {}));
// Error class for SGP4 errors
class Sgp4Error extends Error {
    code;
    constructor(code, message) {
        super(message ?? Sgp4Error.getDefaultMessage(code));
        this.name = 'Sgp4Error';
        this.code = code;
    }
    static getDefaultMessage(code) {
        switch (code) {
            case Sgp4ErrorCode.NO_ERROR:
                return 'No error';
            case Sgp4ErrorCode.MEAN_ELEMENTS_INVALID:
                return 'Mean elements invalid: eccentricity out of range or semi-major axis too small';
            case Sgp4ErrorCode.MEAN_MOTION_NEGATIVE:
                return 'Mean motion is negative';
            case Sgp4ErrorCode.PERT_ELEMENTS_INVALID:
                return 'Perturbed elements invalid: eccentricity out of range';
            case Sgp4ErrorCode.SEMI_LATUS_RECTUM_NEGATIVE:
                return 'Semi-latus rectum is negative';
            case Sgp4ErrorCode.EPOCH_ELEMENTS_SUBORBITAL:
                return 'Epoch elements are sub-orbital';
            case Sgp4ErrorCode.SATELLITE_DECAYED:
                return 'Satellite has decayed';
            default:
                return `Unknown error code: ${code}`;
        }
    }
}


}),
"./src/engine/ootk/src/sgp4/sgp4.ts": 
/*!******************************************!*\
  !*** ./src/engine/ootk/src/sgp4/sgp4.ts ***!
  \******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Sgp4: () => (Sgp4),
  Sgp4GravConstants: () => (Sgp4GravConstants),
  Sgp4Methods: () => (Sgp4Methods)
});
/* ESM import */var _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/Sgp4OpsMode.js */ "./src/engine/ootk/src/enums/Sgp4OpsMode.ts");
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./sgp4-error.js */ "./src/engine/ootk/src/sgp4/sgp4-error.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * This class was ported from the python-sgp4 library by Brandon Rhodes. That library
 * is licensed under the MIT license and he maintains the copyright for that work.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
// NOTE: This file is meant to maintain as much of the original format as possible.
/* eslint-disable complexity */
/* eslint-disable max-statements */
/* eslint-disable max-lines-per-function */
/* eslint-disable max-lines */
/* eslint-disable @typescript-eslint/no-loss-of-precision */




var Sgp4GravConstants;
(function (Sgp4GravConstants) {
    Sgp4GravConstants["wgs72old"] = "wgs72old";
    Sgp4GravConstants["wgs72"] = "wgs72";
    Sgp4GravConstants["wgs84"] = "wgs84";
})(Sgp4GravConstants || (Sgp4GravConstants = {}));
var Sgp4Methods;
(function (Sgp4Methods) {
    Sgp4Methods["NEAR_EARTH"] = "n";
    Sgp4Methods["DEEP_SPACE"] = "d";
})(Sgp4Methods || (Sgp4Methods = {}));
/** Ootk -- Some variables imported from outside the class at the top */
const fasx2 = 0.13130908;
const fasx4 = 2.8843198;
const fasx6 = 0.37448087;
const g22 = 5.7686396;
const g32 = 0.95240898;
const g44 = 1.8014998;
const g52 = 1.050833;
const g54 = 4.4108898;
const rptim = 4.37526908801129966e-3; // Equates to 7.29211514668855e-5 rad/sec
const stepp = 720.0;
const stepn = -720.0;
const step2 = 259200.0;
/*
 *     ----------------------------------------------------------------
 *
 *                               sgp4unit.cpp
 *
 *    this file contains the sgp4 procedures for analytical propagation
 *    of a satellite. the code was originally released in the 1980 and 1986
 *    spacetrack papers. a detailed discussion of the theory and history
 *    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
 *    and kelso.
 *
 *                            companion code for
 *               fundamentals of astrodynamics and applications
 *                                    2013
 *                              by david vallado
 *
 *     (w) 719-573-2600, email dvallado@agi.com, davallado@gmail.com
 *
 *    current :
 *              12 mar 20  david vallado
 *                           chg satnum to string for alpha 5 or 9-digit
 *    changes :
 *               7 dec 15  david vallado
 *                           fix jd, jdfrac
 *               3 nov 14  david vallado
 *                           update to msvs2013 c++
 *              30 aug 10  david vallado
 *                           delete unused variables in initl
 *                           replace pow integer 2, 3 with multiplies for speed
 *               3 nov 08  david vallado
 *                           put returns in for error codes
 *              29 sep 08  david vallado
 *                           fix atime for faster operation in dspace
 *                           add operationmode for afspc (a) or improved (i)
 *                           performance mode
 *              16 jun 08  david vallado
 *                           update small eccentricity check
 *              16 nov 07  david vallado
 *                           misc fixes for better compliance
 *              20 apr 07  david vallado
 *                           misc fixes for constants
 *              11 aug 06  david vallado
 *                           chg lyddane choice back to strn3, constants, misc doc
 *              15 dec 05  david vallado
 *                           misc fixes
 *              26 jul 05  david vallado
 *                           fixes for paper
 *                           note that each fix is preceded by a
 *                           comment with "sgp4fix" and an explanation of
 *                           what was changed
 *              10 aug 04  david vallado
 *                           2nd printing baseline working
 *              14 may 01  david vallado
 *                           2nd edition baseline
 *                     80  norad
 *                           original baseline
 *       ----------------------------------------------------------------
 */
class Sgp4 {
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure angle_SGP4
     *
     *  this procedure calculates the angle between two vectors.  the output is
     *    set to 999999.1 to indicate an undefined value.  be sure to check for
     *    this at the output phase.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    vec1        - vector number 1
     *    vec2        - vector number 2
     *
     *  outputs       :
     *    theta       - angle between the two vectors  -pi to pi
     *
     *  locals        :
     *    temp        - temporary real variable
     *
     *  coupling      :
     *    dot           dot product of two vectors
     * ---------------------------------------------------------------------------
     */
    static angle_(vec1, vec2) {
        const small = 0.00000001;
        const unknown = 999999.1; /** Ootk -- original 'undefined' is protected in JS */
        const magv1 = Sgp4.mag_(vec1);
        const magv2 = Sgp4.mag_(vec2);
        const magnitudeProduct = magv1 * magv2;
        if (magnitudeProduct > small * small) {
            let temp = Sgp4.dot_(vec1, vec2) / (magnitudeProduct);
            // Clamp to [-1, 1] to avoid NaN from floating point errors
            if (temp > 1.0) {
                temp = 1.0;
            }
            if (temp < -1.0) {
                temp = -1.0;
            }
            return Math.acos(temp);
        }
        return unknown;
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function asinh_SGP4
     *
     *  this function evaluates the inverse hyperbolic sine function.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    xval        - angle value                                  any real
     *
     *  outputs       :
     *    arcsinh     - result                                       any real
     *
     *  locals        :
     *    none.
     *
     *  coupling      :
     *    none.
     * ---------------------------------------------------------------------------
     */
    static asinh_(xval) {
        return Math.log(xval + Math.sqrt(xval * xval + 1.0));
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function twoline2rv
     *
     *  this function converts the two line element set character string data to
     *    variables and initializes the sgp4 variables. several intermediate varaibles
     *    and quantities are determined. note that the result is a structure so multiple
     *    satellites can be processed simultaneously without having to reinitialize. the
     *    verification mode is an important option that permits quick checks of any
     *    changes to the underlying technical theory. this option works using a
     *    modified tle file in which the start, stop, and delta time values are
     *    included at the end of the second line of data. this only works with the
     *    verification mode. the catalog mode simply propagates from -1440 to 1440 min
     *    from epoch and is useful when performing entire catalog runs.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs        :
     *    longstr1    - first line of the tle
     *    longstr2    - second line of the tle
     *    typerun     - type of run                    verification 'v', catalog 'c',
     *                                                 manual 'm'
     *    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
     *    opsmode     - mode of operation afspc or improved 'a', 'i'
     *    whichconst  - which set of constants to use  72, 84
     *
     *  outputs       :
     *    satrec      - structure containing all the sgp4 satellite information
     *
     *  coupling      :
     *    getgravconst-
     *    days2mdhms  - conversion of days to month, day, hour, minute, second
     *    jday        - convert day month year hour minute second into julian date
     *    sgp4init    - initialize the sgp4 variables
     *
     *  references    :
     *    norad spacetrack report #3
     *    vallado, crawford, hujsak, kelso  2006
     * ---------------------------------------------------------------------------
     */
    static createSatrec(tleLine1, tleLine2, whichconst = Sgp4GravConstants.wgs72, opsmode = _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.IMPROVED) {
        let year = 0;
        const satrec = {
            error: _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR,
        };
        /*
         * Sgp4fix no longer needed
         * getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
         */
        const xpdotp = 1440.0 / (2.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI); // 229.1831180523293;
        satrec.satnum = tleLine1.substring(2, 7);
        satrec.epochyr = parseInt(tleLine1.substring(18, 20));
        satrec.epochdays = parseFloat(tleLine1.substring(20, 32));
        satrec.ndot = parseFloat(tleLine1.substring(33, 43));
        satrec.nddot = parseFloat(`${tleLine1.substring(44, 45)}.${tleLine1.substring(45, 50)}E${tleLine1.substring(50, 52)}`);
        satrec.bstar = parseFloat(`${tleLine1.substring(53, 54)}.${tleLine1.substring(54, 59)}E${tleLine1.substring(59, 61)}`);
        satrec.inclo = parseFloat(tleLine2.substring(8, 16));
        satrec.nodeo = parseFloat(tleLine2.substring(17, 25));
        satrec.ecco = parseFloat(`.${tleLine2.substring(26, 33)}`);
        satrec.argpo = parseFloat(tleLine2.substring(34, 42));
        satrec.mo = parseFloat(tleLine2.substring(43, 51));
        satrec.no = parseFloat(tleLine2.substring(52, 63));
        // ---- find no, ndot, nddot ----
        satrec.no /= xpdotp; //   Rad/min
        /** Ootk -- nexp and ibexp are calculated above using template literals */
        /*
         * Satrec.nddot = satrec.nddot * Math.pow(10.0, nexp);
         * satrec.bstar = satrec.bstar * Math.pow(10.0, ibexp);
         */
        /*
         * ---- convert to sgp4 units ----
         * satrec.a = (satrec.no * tumin) ** (-2.0 / 3.0);
         */
        /** Ootk -- Not sure why the following two lines are added. 1st and 2nd derivatives aren't even used anymore */
        /*
         * Satrec.ndot /= xpdotp * 1440.0; // ? * minperday
         * satrec.nddot /= xpdotp * 1440.0 * 1440;
         */
        // ---- find standard orbital elements ----
        satrec.inclo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.nodeo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.argpo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.mo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        /*
         * Sgp4fix not needed here
         * satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
         * satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
         */
        /*
         * ----------------------------------------------------------------
         * find sgp4epoch time of element set
         * remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
         * and minutes from the epoch (time)
         * ----------------------------------------------------------------
         */
        /*
         * ---------------- temp fix for years from 1957-2056 -------------------
         * --------- correct fix will occur when year is 4-digit in tle ---------
         */
        if (satrec.epochyr < 57) {
            year = satrec.epochyr + 2000;
        }
        else {
            year = satrec.epochyr + 1900;
        }
        const { mon, day, hr, min, sec } = Sgp4.days2mdhms(year, satrec.epochdays);
        const jdayRes = Sgp4.jday(year, mon, day, hr, min, sec);
        satrec.jdsatepoch = jdayRes.jd + jdayRes.jdFrac;
        //  ---------------- initialize the orbit at sgp4epoch -------------------
        Sgp4.sgp4init_(satrec, {
            whichconst,
            opsmode,
            satn: satrec.satnum,
            epoch: satrec.jdsatepoch - 2433281.5,
            xbstar: satrec.bstar,
            xecco: satrec.ecco,
            xargpo: satrec.argpo,
            xinclo: satrec.inclo,
            xndot: satrec.ndot,
            xnddot: satrec.nddot,
            xmo: satrec.mo,
            xno: satrec.no,
            xnodeo: satrec.nodeo,
        });
        return satrec;
    }
    static createSatrecFromOmm(omm, whichconst = Sgp4GravConstants.wgs72, opsmode = _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.IMPROVED) {
        let year = 0;
        const satrec = {
            error: _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR,
        };
        const xpdotp = 1440.0 / (2.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI); // 229.1831180523293;
        satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR;
        satrec.satnum = omm.NORAD_CAT_ID;
        satrec.epochyr = parseInt(omm.epoch.year.toString().slice(-2));
        satrec.epochdays = omm.epoch.doy;
        satrec.ndot = parseFloat(omm.MEAN_MOTION_DOT);
        satrec.nddot = parseFloat(omm.MEAN_MOTION_DDOT);
        satrec.bstar = parseFloat(omm.BSTAR);
        satrec.inclo = parseFloat(omm.INCLINATION);
        satrec.nodeo = parseFloat(omm.RA_OF_ASC_NODE);
        satrec.ecco = parseFloat(omm.ECCENTRICITY);
        satrec.argpo = parseFloat(omm.ARG_OF_PERICENTER);
        satrec.mo = parseFloat(omm.MEAN_ANOMALY);
        satrec.no = parseFloat(omm.MEAN_MOTION);
        // ---- find no, ndot, nddot ----
        satrec.no /= xpdotp; //   Rad/min
        /** Ootk -- nexp and ibexp are calculated above using template literals */
        /*
         * Satrec.nddot = satrec.nddot * Math.pow(10.0, nexp);
         * satrec.bstar = satrec.bstar * Math.pow(10.0, ibexp);
         */
        /*
         * ---- convert to sgp4 units ----
         * satrec.a = (satrec.no * tumin) ** (-2.0 / 3.0);
         */
        /** Ootk -- Not sure why the following two lines are added. 1st and 2nd derivatives aren't even used anymore */
        /*
         * Satrec.ndot /= xpdotp * 1440.0; // ? * minperday
         * satrec.nddot /= xpdotp * 1440.0 * 1440;
         */
        // ---- find standard orbital elements ----
        satrec.inclo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.nodeo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.argpo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        satrec.mo *= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD;
        /*
         * Sgp4fix not needed here
         * satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
         * satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
         */
        /*
         * ----------------------------------------------------------------
         * find sgp4epoch time of element set
         * remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
         * and minutes from the epoch (time)
         * ----------------------------------------------------------------
         */
        /*
         * ---------------- temp fix for years from 1957-2056 -------------------
         * --------- correct fix will occur when year is 4-digit in tle ---------
         */
        if (satrec.epochyr < 57) {
            year = satrec.epochyr + 2000;
        }
        else {
            year = satrec.epochyr + 1900;
        }
        const { mon, day, hr, min, sec } = Sgp4.days2mdhms(year, satrec.epochdays);
        const jdayRes = Sgp4.jday(year, mon, day, hr, min, sec);
        satrec.jdsatepoch = jdayRes.jd + jdayRes.jdFrac;
        //  ---------------- initialize the orbit at sgp4epoch -------------------
        Sgp4.sgp4init_(satrec, {
            whichconst,
            opsmode,
            satn: satrec.satnum,
            epoch: satrec.jdsatepoch - 2433281.5,
            xbstar: satrec.bstar,
            xecco: satrec.ecco,
            xargpo: satrec.argpo,
            xinclo: satrec.inclo,
            xndot: satrec.ndot,
            xnddot: satrec.nddot,
            xmo: satrec.mo,
            xno: satrec.no,
            xnodeo: satrec.nodeo,
        });
        return satrec;
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure cross_SGP4
     *
     *  this procedure crosses two vectors.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    vec1        - vector number 1
     *    vec2        - vector number 2
     *
     *  outputs       :
     *    outvec      - vector result of a x b
     *
     *  locals        :
     *    none.
     *
     *  coupling      :
     *    mag           magnitude of a vector
     * ----------------------------------------------------------------------------
     */
    static cross_(vec1, vec2) {
        return [
            vec1[1] * vec2[2] - vec1[2] * vec2[1],
            vec1[2] * vec2[0] - vec1[0] * vec2[2],
            vec1[0] * vec2[1] - vec1[1] * vec2[0],
        ];
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure days2mdhms
     *
     *  this procedure converts the day of the year, days, to the equivalent month
     *    day, hour, minute and second.
     *
     *  algorithm     : set up array for the number of days per month
     *                  find leap year - use 1900 because 2000 is a leap year
     *                  loop through a temp value while the value is < the days
     *                  perform int conversions to the correct day and month
     *                  convert remainder into h m s using type conversions
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    year        - year                           1900 .. 2100
     *    days        - julian day of the year         0.0  .. 366.0
     *
     *  outputs       :
     *    mon         - month                          1 .. 12
     *    day         - day                            1 .. 28,29,30,31
     *    hr          - hour                           0 .. 23
     *    min         - minute                         0 .. 59
     *    sec         - second                         0.0 .. 59.999
     *
     *  locals        :
     *    dayofyr     - day of year
     *    temp        - temporary extended values
     *    inttemp     - temporary int value
     *    i           - index
     *    lmonth[13]  - int array containing the number of days per month
     *
     *  coupling      :
     *    none.
     * ---------------------------------------------------------------------------
     */
    static days2mdhms(year, days) {
        const lmonth = [31, year % 4 === 0 ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
        const dayofyr = Math.floor(days);
        //  ----------------- find month and day of month ----------------
        /** Ootk -- Incorporated in the above declaration */
        /*
         * If ((year % 4) == 0)
         * lmonth[2] = 29;
         */
        let i = 1;
        let inttemp = 0;
        while (dayofyr > inttemp + (lmonth[i - 1]) && i < 12) {
            inttemp += (lmonth[i - 1]);
            i += 1;
        }
        const mon = i;
        const day = dayofyr - inttemp;
        //  ----------------- find hours minutes and seconds -------------
        let temp = (days - dayofyr) * 24.0;
        const hr = Math.floor(temp);
        temp = (temp - hr) * 60.0;
        const min = Math.floor(temp);
        const sec = (temp - min) * 60.0;
        return {
            mon,
            day,
            hr,
            min,
            sec,
        };
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function dot_SGP4
     *
     *  this function finds the dot product of two vectors.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    vec1        - vector number 1
     *    vec2        - vector number 2
     *
     *  outputs       :
     *    dot         - result
     *
     *  locals        :
     *    none.
     *
     *  coupling      :
     *    none.
     * ---------------------------------------------------------------------------
     */
    static dot_(v1, v2) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function gstime
     *
     *  this function finds the greenwich sidereal time.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    jdut1       - julian date in ut1             days from 4713 bc
     *
     *  outputs       :
     *    gstime      - greenwich sidereal time        0 to 2PI rad
     *
     *  locals        :
     *    temp        - temporary variable for doubles   rad
     *    tut1        - julian centuries from the
     *                  jan 1, 2000 12 h epoch (ut1)
     *
     *  coupling      :
     *    none
     *
     *  references    :
     *    vallado       2004, 191, eq 3-45
     * ---------------------------------------------------------------------------
     */
    static gstime(jdut1) {
        const tut1 = (jdut1 - 2451545.0) / 36525.0;
        let temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
            (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841; // Sec
        temp = ((temp * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD) / 240.0) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU; // 360/86400 = 1/240, to deg, to rad
        //  ------------------------ check quadrants ---------------------
        if (temp < 0.0) {
            temp += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        }
        return temp;
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure invjday
     *
     *  this procedure finds the year, month, day, hour, minute and second
     *  given the julian date. tu can be ut1, tdt, tdb, etc.
     *
     *  algorithm     : set up starting values
     *                  find leap year - use 1900 because 2000 is a leap year
     *                  find the elapsed days through the year in a loop
     *                  call routine to find each individual value
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    jd          - julian date                    days from 4713 bc
     *    jdfrac      - julian date fraction into day  days from 4713 bc
     *
     *  outputs       :
     *    year        - year                           1900 .. 2100
     *    mon         - month                          1 .. 12
     *    day         - day                            1 .. 28,29,30,31
     *    hr          - hour                           0 .. 23
     *    min         - minute                         0 .. 59
     *    sec         - second                         0.0 .. 59.999
     *
     *  locals        :
     *    days        - day of year plus fractional
     *                  portion of a day               days
     *    tu          - julian centuries from 0 h
     *                  jan 0, 1900
     *    temp        - temporary double values
     *    leapyrs     - number of leap years from 1900
     *
     *  coupling      :
     *    days2mdhms  - finds month, day, hour, minute and second given days and year
     *
     *  references    :
     *    vallado       2013, 203, alg 22, ex 3-13
     * ---------------------------------------------------------------------------
     */
    static invjday(jd, jdfrac) {
        let leapyrs;
        let days;
        // Check jdfrac for multiple days
        if (Math.abs(jdfrac) >= 1.0) {
            jd += Math.floor(jdfrac);
            jdfrac -= Math.floor(jdfrac);
        }
        // Check for fraction of a day included in the jd
        const dt = jd - Math.floor(jd) - 0.5;
        if (Math.abs(dt) > 0.00000001) {
            jd -= dt;
            jdfrac += dt;
        }
        /* --------------- find year and days of the year --------------- */
        const temp = jd - 2415019.5;
        const tu = temp / 365.25;
        let year = 1900 + Math.floor(tu);
        leapyrs = Math.floor((year - 1901) * 0.25);
        days = Math.floor(temp - ((year - 1900) * 365.0 + leapyrs));
        /* ------------ check for case of beginning of a year ----------- */
        if (days + jdfrac < 1.0) {
            year -= 1;
            leapyrs = Math.floor((year - 1901) * 0.25);
            days = Math.floor(temp - ((year - 1900) * 365.0 + leapyrs));
        }
        /* ----------------- find remaining data  ------------------------- */
        const { mon, day, hr, min, sec } = Sgp4.days2mdhms(year, days + jdfrac);
        return {
            year,
            mon,
            day,
            hr,
            min,
            sec,
        };
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure jday
     *
     *  this procedure finds the julian date given the year, month, day, and time.
     *    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
     *
     *  algorithm     : calculate the answer in one step for efficiency
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    year        - year                           1900 .. 2100
     *    mon         - month                          1 .. 12
     *    day         - day                            1 .. 28,29,30,31
     *    hr          - universal time hour            0 .. 23
     *    min         - universal time min             0 .. 59
     *    sec         - universal time sec             0.0 .. 59.999
     *
     *  outputs       :
     *    jd          - julian date                    days from 4713 bc
     *    jdfrac      - julian date fraction into day  days from 4713 bc
     *
     *  locals        :
     *    none.
     *
     *  coupling      :
     *    none.
     *
     *  references    :
     *    vallado       2013, 183, alg 14, ex 3-4
     *
     * ---------------------------------------------------------------------------
     */
    static jday(year, mon = 0, day = 0, hr = 0, min = 0, sec = 0, ms = 0) {
        if (year instanceof Date) {
            mon = year.getUTCMonth() + 1;
            day = year.getUTCDate();
            hr = year.getUTCHours();
            min = year.getUTCMinutes();
            sec = year.getUTCSeconds();
            ms = year.getUTCMilliseconds();
            year = year.getUTCFullYear();
        }
        let jd = 367.0 * year -
            Math.floor(7 * (year + Math.floor((mon + 9) / 12.0)) * 0.25) +
            Math.floor((275 * mon) / 9.0) +
            day +
            1721013.5; // Use - 678987.0 to go to mjd directly
        let jdFrac = (ms / 1000 + sec + min * 60.0 + hr * 3600.0) / 86400.0;
        // Check that the day and fractional day are correct
        if (Math.abs(jdFrac) > 1.0) {
            const dtt = Math.floor(jdFrac);
            jd += dtt;
            jdFrac -= dtt;
        }
        // - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
        return { jd, jdFrac };
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function mag
     *
     *  this procedure finds the magnitude of a vector.
     *
     *  author        : david vallado                  719-573-2600    1 mar 2001
     *
     *  inputs          description                    range / units
     *    vec         - vector
     *
     *  outputs       :
     *    mag         - answer
     *
     *  locals        :
     *    none.
     *
     *  coupling      :
     *    none.
     * ---------------------------------------------------------------------------
     */
    static mag_(v) {
        return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function newtonnu_SGP4
     *
     *  this function solves keplers equation when the true anomaly is known.
     *    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
     *    the parabolic limit at 168ø is arbitrary. the hyperbolic anomaly is also
     *    limited. the hyperbolic sine is used because it's not double valued.
     *
     *  author        : david vallado                  719-573-2600   27 may 2002
     *
     *  revisions
     *    vallado     - fix small                                     24 sep 2002
     *
     *  inputs          description                    range / units
     *    ecc         - eccentricity                   0.0  to
     *    nu          - true anomaly                   -2pi to 2pi rad
     *
     *  outputs       :
     *    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 ø
     *    m           - mean anomaly                   0.0  to 2pi rad       151.7425 ø
     *
     *  locals        :
     *    e1          - eccentric anomaly, next value  rad
     *    sine        - sine of e
     *    cose        - cosine of e
     *    ktr         - index
     *
     *  coupling      :
     *    asinh       - arc hyperbolic sine
     *
     *  references    :
     *    vallado       2013, 77, alg 5
     * ---------------------------------------------------------------------------
     */
    static newtonnu_(ecc, nu) {
        // ---------------------  implementation   ---------------------
        let e0 = 999999.9;
        let m = 999999.9;
        const small = 0.00000001;
        if (Math.abs(ecc) < small) {
            // --------------------------- circular ------------------------
            m = nu;
            e0 = nu;
        }
        else if (ecc < 1.0 - small) {
            // ---------------------- elliptical -----------------------
            const sine = (Math.sqrt(1.0 - ecc * ecc) * Math.sin(nu)) / (1.0 + ecc * Math.cos(nu));
            const cose = (ecc + Math.cos(nu)) / (1.0 + ecc * Math.cos(nu));
            e0 = Math.atan2(sine, cose);
            m = e0 - ecc * Math.sin(e0);
        }
        else if (ecc > 1.0 + small) {
            // -------------------- hyperbolic  --------------------
            if (ecc > 1.0 && Math.abs(nu) + 0.00001 < _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI - Math.acos(1.0 / ecc)) {
                const sine = (Math.sqrt(ecc * ecc - 1.0) * Math.sin(nu)) / (1.0 + ecc * Math.cos(nu));
                e0 = Sgp4.asinh_(sine);
                m = ecc * Sgp4.sinh_(e0) - e0;
            }
        }
        else if (Math.abs(nu) < (168.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI) / 180.0) {
            // ----------------- parabolic ---------------------
            e0 = Math.tan(nu * 0.5);
            m = e0 + (e0 * e0 * e0) / 3.0;
        }
        if (ecc < 1.0) {
            m %= 2.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
            if (m < 0.0) {
                m += 2.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
            }
            e0 %= 2.0 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
        }
        return {
            e0,
            m,
        };
    }
    /*
     *----------------------------------------------------------------------------
     *
     *                             procedure sgp4
     *
     *  this procedure is the sgp4 prediction model from space command. this is an
     *    updated and combined version of sgp4 and sdp4, which were originally
     *    published separately in spacetrack report //3. this version follows the
     *    methodology from the aiaa paper (2006) describing the history and
     *    development of the code.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    satrec  - initialised structure from sgp4init() call.
     *    tsince  - time since epoch (minutes)
     *
     *  outputs       :
     *    r           - position vector                     km
     *    v           - velocity                            km/sec
     *  return code - non-zero on error.
     *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
     *                   2 - mean motion less than 0.0
     *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
     *                   4 - semi-latus rectum < 0.0
     *                   5 - epoch elements are sub-orbital
     *                   6 - satellite has decayed
     *
     *  locals        :
     *    am          -
     *    axnl, aynl        -
     *    betal       -
     *    cosim   , sinim   , cosomm  , sinomm  , cnod    , snod    , cos2u   ,
     *    sin2u   , coseo1  , sineo1  , cosi    , sini    , cosip   , sinip   ,
     *    cosisq  , cossu   , sinsu   , cosu    , sinu
     *    delm        -
     *    delomg      -
     *    dndt        -
     *    eccm        -
     *    emsq        -
     *    ecose       -
     *    el2         -
     *    eo1         -
     *    eccp        -
     *    esine       -
     *    argpm       -
     *    argpp       -
     *    omgadf      -
     *    pl          -
     *    r           -
     *    rtemsq      -
     *    rdotl       -
     *    rl          -
     *    rvdot       -
     *    rvdotl      -
     *    su          -
     *    t2  , t3   , t4    , tc
     *    tem5, temp , temp1 , temp2  , tempa  , tempe  , templ
     *    u   , ux   , uy    , uz     , vx     , vy     , vz
     *    inclm       - inclination
     *    mm          - mean anomaly
     *    nm          - mean motion
     *    nodem       - right asc of ascending node
     *    xinc        -
     *    xincp       -
     *    xl          -
     *    xlm         -
     *    mp          -
     *    xmdf        -
     *    xmx         -
     *    xmy         -
     *    nodedf      -
     *    xnode       -
     *    nodep       -
     *    np          -
     *
     *  coupling      :
     *    getgravconst-
     *    dpper
     *    dspace
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report //3 1980
     *    hoots, norad spacetrack report //6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static propagate(satrec, tsince) {
        /* ------------------ set mathematical constants --------------- */
        /*
         * Sgp4fix divisor for divide by zero check on inclination
         * the old check used 1.0 + cos(PI-1.0e-9), but then compared it to
         * 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
         */
        /*
         * Sgp4fix identify constants and allow alternate values
         * getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
         */
        const { xke, j2, j3oj2, vkmpersec } = satrec;
        // --------------------- clear sgp4 error flag -----------------
        satrec.t = tsince;
        satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR;
        //  ------- update for secular gravity and atmospheric drag -----
        const xmdf = satrec.mo + satrec.mdot * satrec.t;
        const argpdf = satrec.argpo + satrec.argpdot * satrec.t;
        const nodedf = satrec.nodeo + satrec.nodedot * satrec.t;
        let argpm = argpdf;
        let mm = xmdf;
        const t2 = satrec.t * satrec.t;
        let nodem = nodedf + satrec.nodecf * t2;
        let tempa = 1.0 - satrec.cc1 * satrec.t;
        let tempe = satrec.bstar * satrec.cc4 * satrec.t;
        let templ = satrec.t2cof * t2;
        if (!satrec.isimp) {
            const delomg = satrec.omgcof * satrec.t;
            //  Sgp4fix use mutliply for speed instead of pow
            const delmtemp = 1.0 + satrec.eta * Math.cos(xmdf);
            const delm = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo);
            const temp = delomg + delm;
            mm = xmdf + temp;
            argpm = argpdf - temp;
            const t3 = t2 * satrec.t;
            const t4 = t3 * satrec.t;
            tempa = tempa - satrec.d2 * t2 - satrec.d3 * t3 - satrec.d4 * t4;
            tempe += satrec.bstar * satrec.cc5 * (Math.sin(mm) - satrec.sinmao);
            templ = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof + satrec.t * satrec.t5cof);
        }
        let nm = satrec.no;
        let em = satrec.ecco;
        let inclm = satrec.inclo;
        if (satrec.method === Sgp4Methods.DEEP_SPACE) {
            [em, argpm, inclm, mm, nodem, nm] = Sgp4.dspace_(em, argpm, inclm, mm, nodem, nm, satrec);
        }
        if (nm <= 0.0) {
            satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.MEAN_MOTION_NEGATIVE;
            return { position: false, velocity: false };
        }
        const am = (xke / nm) ** _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.x2o3 * tempa * tempa;
        nm = xke / am ** 1.5;
        em -= tempe;
        /*
         * Fix tolerance for error recognition
         * sgp4fix am is fixed from the previous nm check
         */
        /* istanbul ignore next | This is no longer possible*/
        if (em >= 1.0 || em < -0.001) {
            satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.MEAN_MOTION_NEGATIVE;
            return { position: false, velocity: false };
        }
        //  Sgp4fix fix tolerance to avoid a divide by zero
        if (em < 1.0e-6) {
            em = 1.0e-6;
        }
        mm += satrec.no * templ;
        let xlm = mm + argpm + nodem;
        nodem %= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        argpm %= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        xlm %= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        mm = (xlm - argpm - nodem) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        /*
         * Sgp4fix recover singly averaged mean elements
         * satrec.am = am;
         * satrec.em = em;
         * satrec.im = inclm;
         * satrec.Om = nodem;
         * satrec.om = argpm;
         * satrec.mm = mm;
         * satrec.nm = nm;
         */
        // ----------------- compute extra mean quantities -------------
        const sinim = Math.sin(inclm);
        const cosim = Math.cos(inclm);
        // -------------------- add lunar-solar periodics --------------
        let ep = em;
        let xincp = inclm;
        let argpp = argpm;
        let nodep = nodem;
        let mp = mm;
        let sinip = sinim;
        let cosip = cosim;
        if (satrec.method === Sgp4Methods.DEEP_SPACE) {
            const dpperParameters = {
                inclo: satrec.inclo,
                init: false,
                ep,
                inclp: xincp,
                nodep,
                argpp,
                mp,
                opsmode: satrec.operationmode,
                satrec,
            };
            ({ ep, nodep, argpp, mp, inclp: xincp } = Sgp4.dpper_(dpperParameters));
            if (xincp < 0.0) {
                xincp = -xincp;
                nodep += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
                argpp -= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
            }
            if (ep < 0.0 || ep > 1.0) {
                satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.PERT_ELEMENTS_INVALID;
                return { position: false, velocity: false };
            }
        }
        //  -------------------- long period periodics ------------------
        if (satrec.method === Sgp4Methods.DEEP_SPACE) {
            sinip = Math.sin(xincp);
            cosip = Math.cos(xincp);
            satrec.aycof = -0.5 * j3oj2 * sinip;
            //  Sgp4fix for divide by zero for xincp = 180 deg
            if (Math.abs(cosip + 1.0) > 1.5e-12) {
                satrec.xlcof = (-0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip)) / (1.0 + cosip);
            }
            else {
                satrec.xlcof = (-0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip)) / _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.temp4;
            }
        }
        const axnl = ep * Math.cos(argpp);
        let temp = 1.0 / (am * (1.0 - ep * ep));
        const aynl = ep * Math.sin(argpp) + temp * satrec.aycof;
        const xl = mp + argpp + nodep + temp * satrec.xlcof * axnl;
        // --------------------- solve kepler's equation ---------------
        const u = (xl - nodep) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        let eo1 = u;
        let tem5 = 9999.9;
        let ktr = 1;
        /*
         *    Sgp4fix for kepler iteration
         *    the following iteration needs better limits on corrections
         */
        let coseo1 = 0;
        let sineo1 = 0;
        while (Math.abs(tem5) >= 1.0e-12 && ktr <= 10) {
            sineo1 = Math.sin(eo1);
            coseo1 = Math.cos(eo1);
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if (Math.abs(tem5) >= 0.95) {
                if (tem5 > 0.0) {
                    tem5 = 0.95;
                }
                else {
                    tem5 = -0.95;
                }
            }
            eo1 += tem5;
            ktr += 1;
        }
        //  ------------- short period preliminary quantities -----------
        const ecose = axnl * coseo1 + aynl * sineo1;
        const esine = axnl * sineo1 - aynl * coseo1;
        const el2 = axnl * axnl + aynl * aynl;
        const pl = am * (1.0 - el2);
        if (pl < 0.0) {
            satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.SEMI_LATUS_RECTUM_NEGATIVE;
            return { position: false, velocity: false };
        }
        const rl = am * (1.0 - ecose);
        const rdotl = (Math.sqrt(am) * esine) / rl;
        const rvdotl = Math.sqrt(pl) / rl;
        const betal = Math.sqrt(1.0 - el2);
        temp = esine / (1.0 + betal);
        const sinu = (am / rl) * (sineo1 - aynl - axnl * temp);
        const cosu = (am / rl) * (coseo1 - axnl + aynl * temp);
        let su = Math.atan2(sinu, cosu);
        const sin2u = (cosu + cosu) * sinu;
        const cos2u = 1.0 - 2.0 * sinu * sinu;
        temp = 1.0 / pl;
        const temp1 = 0.5 * j2 * temp;
        const temp2 = temp1 * temp;
        // -------------- update for short period periodics ------------
        if (satrec.method === Sgp4Methods.DEEP_SPACE) {
            const cosisq = cosip * cosip;
            satrec.con41 = 3.0 * cosisq - 1.0;
            satrec.x1mth2 = 1.0 - cosisq;
            satrec.x7thm1 = 7.0 * cosisq - 1.0;
        }
        const mrt = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + 0.5 * temp1 * satrec.x1mth2 * cos2u;
        /** Moved this up to reduce unnecessary computation if you are going to return false anyway */
        // Sgp4fix for decaying satellites
        if (mrt < 1.0) {
            satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.SATELLITE_DECAYED;
            return {
                position: false,
                velocity: false,
            };
        }
        su -= 0.25 * temp2 * satrec.x7thm1 * sin2u;
        const xnode = nodep + 1.5 * temp2 * cosip * sin2u;
        const xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
        const mvt = rdotl - (nm * temp1 * satrec.x1mth2 * sin2u) / xke;
        const rvdot = rvdotl + (nm * temp1 * (satrec.x1mth2 * cos2u + 1.5 * satrec.con41)) / xke;
        // --------------------- orientation vectors -------------------
        const sinsu = Math.sin(su);
        const cossu = Math.cos(su);
        const snod = Math.sin(xnode);
        const cnod = Math.cos(xnode);
        const sini = Math.sin(xinc);
        const cosi = Math.cos(xinc);
        const xmx = -snod * cosi;
        const xmy = cnod * cosi;
        const ux = xmx * sinsu + cnod * cossu;
        const uy = xmy * sinsu + snod * cossu;
        const uz = sini * sinsu;
        const vx = xmx * cossu - cnod * sinsu;
        const vy = xmy * cossu - snod * sinsu;
        const vz = sini * cossu;
        // --------- position and velocity (in km and km/sec) ----------
        const r = {
            x: mrt * ux * satrec.radiusearthkm,
            y: mrt * uy * satrec.radiusearthkm,
            z: mrt * uz * satrec.radiusearthkm,
        };
        const v = {
            x: (mvt * ux + rvdot * vx) * vkmpersec,
            y: (mvt * uy + rvdot * vy) * vkmpersec,
            z: (mvt * uz + rvdot * vz) * vkmpersec,
        };
        return {
            position: r,
            velocity: v,
        };
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function rv2coe_SGP4
     *
     *  this function finds the classical orbital elements given the geocentric
     *    equatorial position and velocity vectors.
     *
     *  author        : david vallado                  719-573-2600   21 jun 2002
     *
     *  revisions
     *    vallado     - fix special cases                              5 sep 2002
     *    vallado     - delete extra check in inclination code        16 oct 2002
     *    vallado     - add constant file use                         29 jun 2003
     *    vallado     - add mu                                         2 apr 2007
     *
     *  inputs          description                    range / units
     *    r           - ijk position vector            km
     *    v           - ijk velocity vector            km / s
     *    mu          - gravitational parameter        km3 / s2
     *
     *  outputs       :
     *    p           - semilatus rectum               km
     *    a           - semimajor axis                 km
     *    ecc         - eccentricity
     *    incl        - inclination                    0.0  to pi rad
     *    omega       - right ascension of ascending node    0.0  to 2pi rad
     *    argp        - argument of perigee            0.0  to 2pi rad
     *    nu          - true anomaly                   0.0  to 2pi rad
     *    m           - mean anomaly                   0.0  to 2pi rad
     *    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
     *    truelon     - true longitude            (ce) 0.0  to 2pi rad
     *    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
     *
     *  locals        :
     *    hbar        - angular momentum h vector      km2 / s
     *    ebar        - eccentricity     e vector
     *    nbar        - line of nodes    n vector
     *    c1          - v**2 - u/r
     *    rdotv       - r dot v
     *    hk          - hk unit vector
     *    sme         - specfic mechanical energy      km2 / s2
     *    i           - index
     *    e           - eccentric, parabolic,
     *                  hyperbolic anomaly             rad
     *    temp        - temporary variable
     *    typeorbit   - type of orbit                  ee, ei, ce, ci
     *
     *  coupling      :
     *    mag         - magnitude of a vector
     *    cross       - cross product of two vectors
     *    angle       - find the angle between two vectors
     *    newtonnu    - find the mean anomaly
     *
     *  references    :
     *    vallado       2013, 113, alg 9, ex 2-5
     * ---------------------------------------------------------------------------
     */
    static rv2coe(r, v, mus) {
        const nbar = [0, 0, 0];
        const ebar = [0, 0, 0];
        let p;
        let a;
        let ecc;
        let incl;
        let omega;
        let argp;
        let nu;
        let m = 0;
        let arglat;
        let truelon;
        let lonper;
        let rdotv;
        let magn;
        let hk;
        let sme;
        let i;
        /*
         *  Switch this to an integer msvs seems to have probelms with this and strncpy_s
         * char typeorbit[2];
         */
        let typeorbit;
        /*
         * Here
         * typeorbit = 1 = 'ei'
         * typeorbit = 2 = 'ce'
         * typeorbit = 3 = 'ci'
         * typeorbit = 4 = 'ee'
         */
        const halfpi = 0.5 * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI;
        const small = 0.00000001;
        const unknown = 999999.1; /** Ootk -- original undefined is illegal in JS */
        const infinite = 999999.9;
        // -------------------------  implementation   -----------------
        const magr = Sgp4.mag_(r);
        const magv = Sgp4.mag_(v);
        // ------------------  find h n and e vectors   ----------------
        const hbar = Sgp4.cross_(r, v);
        const magh = Sgp4.mag_(hbar);
        if (magh > small) {
            nbar[0] = -hbar[1];
            nbar[1] = hbar[0];
            nbar[2] = 0.0;
            magn = Sgp4.mag_(nbar);
            const c1 = magv * magv - mus / magr;
            rdotv = Sgp4.dot_(r, v);
            for (i = 0; i <= 2; i++) {
                ebar[i] = (c1 * (r[i]) - rdotv * (v[i])) / mus;
            }
            ecc = Sgp4.mag_(ebar);
            // ------------  find a e and semi-latus rectum   ----------
            sme = magv * magv * 0.5 - mus / magr;
            if (Math.abs(sme) > small) {
                a = -mus / (2.0 * sme);
            }
            else {
                a = infinite;
            }
            p = (magh * magh) / mus;
            // -----------------  find inclination   -------------------
            hk = hbar[2] / magh;
            incl = Math.acos(hk);
            /*
             * --------  determine type of orbit for later use  --------
             * ------ elliptical, parabolic, hyperbolic inclined -------
             */
            typeorbit = 1;
            if (ecc < small) {
                // ----------------  circular equatorial ---------------
                if (incl < small || Math.abs(incl - _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI) < small) {
                    typeorbit = 2;
                }
                else {
                    // --------------  circular inclined ---------------
                    typeorbit = 3;
                }
                // - elliptical, parabolic, hyperbolic equatorial --
            }
            else if (incl < small || Math.abs(incl - _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI) < small) {
                typeorbit = 4;
            }
            // ----------  find right ascension of the ascending node ------------
            if (magn > small) {
                let temp = nbar[0] / magn;
                if (Math.abs(temp) > 1.0) {
                    temp = Sgp4.sgn_(temp);
                }
                omega = Math.acos(temp);
                if (nbar[1] < 0.0) {
                    omega = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - omega;
                }
            }
            else {
                omega = unknown;
            }
            // ---------------- find argument of perigee ---------------
            if (typeorbit === 1) {
                argp = Sgp4.angle_(nbar, ebar);
                if (ebar[2] < 0.0) {
                    argp = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - argp;
                }
            }
            else {
                argp = unknown;
            }
            // ------------  find true anomaly at epoch    -------------
            if (typeorbit === 1 || typeorbit === 4) {
                nu = Sgp4.angle_(ebar, r);
                if (rdotv < 0.0) {
                    nu = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - nu;
                }
            }
            else {
                nu = unknown;
            }
            // ----  find argument of latitude - circular inclined -----
            if (typeorbit === 3) {
                arglat = Sgp4.angle_(nbar, r);
                if (r[2] < 0.0) {
                    arglat = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - arglat;
                }
                m = arglat;
            }
            else {
                arglat = unknown;
            }
            if (ecc > small && typeorbit === 4) {
                let temp = ebar[0] / ecc;
                if (Math.abs(temp) > 1.0) {
                    temp = Sgp4.sgn_(temp);
                }
                lonper = Math.acos(temp);
                if (ebar[1] < 0.0) {
                    lonper = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - lonper;
                }
                if (incl > halfpi) {
                    lonper = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - lonper;
                }
            }
            else {
                lonper = unknown;
            }
            // -------- find true longitude - circular equatorial ------
            if (magr > small && typeorbit === 2) {
                let temp = r[0] / magr;
                if (Math.abs(temp) > 1.0) {
                    temp = Sgp4.sgn_(temp);
                }
                truelon = Math.acos(temp);
                if (r[1] < 0.0) {
                    truelon = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - truelon;
                }
                if (incl > halfpi) {
                    truelon = _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU - truelon;
                }
                m = truelon;
            }
            else {
                truelon = unknown;
            }
            // ------------ find mean anomaly for all orbits -----------
            if (typeorbit === 1 || typeorbit === 4) {
                m = Sgp4.newtonnu_(ecc, nu).m;
            }
        }
        else {
            p = unknown;
            a = unknown;
            ecc = unknown;
            incl = unknown;
            omega = unknown;
            argp = unknown;
            nu = unknown;
            m = unknown;
            arglat = unknown;
            truelon = unknown;
            lonper = unknown;
        }
        return {
            p,
            a,
            ecc,
            incl,
            omega,
            argp,
            nu,
            m,
            arglat,
            truelon,
            lonper,
        };
    }
    /**
     * Determines the sign of a given number.
     * @param x - The input number to evaluate.
     * @returns `-1.0` if the input number is less than `0.0`, otherwise `1.0`.
     */
    static sgn_(x) {
        if (x < 0.0) {
            return -1.0;
        }
        return 1.0;
    }
    /**
     * Computes the hyperbolic sine of a given number.
     *
     * The hyperbolic sine is calculated using the formula:
     * sinh(x) = (e^x - e^(-x)) / 2
     * @param x - The input number for which to calculate the hyperbolic sine.
     * @returns The hyperbolic sine of the input number.
     */
    static sinh_(x) {
        return (Math.exp(x) - Math.exp(-x)) / 2;
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           procedure dpper
     *
     *  this procedure provides deep space long period periodic contributions
     *    to the mean elements.  by design, these periodics are zero at epoch.
     *    this used to be dscom which included initialization, but it's really a
     *    recurring function.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    e3          -
     *    ee2         -
     *    peo         -
     *    pgho        -
     *    pho         -
     *    PInco       -
     *    plo         -
     *    se2 , se3 , sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4 -
     *    t           -
     *    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
     *    zmol        -
     *    zmos        -
     *    ep          - eccentricity                           0.0 - 1.0
     *    inclo       - inclination - needed for lyddane modification
     *    nodep       - right ascension of ascending node
     *    argpp       - argument of perigee
     *    mp          - mean anomaly
     *
     *  outputs       :
     *    ep          - eccentricity                           0.0 - 1.0
     *    inclp       - inclination
     *    nodep        - right ascension of ascending node
     *    argpp       - argument of perigee
     *    mp          - mean anomaly
     *
     *  locals        :
     *    alfdp       -
     *    betdp       -
     *    cosip  , sinip  , cosop  , sinop  ,
     *    dalf        -
     *    dbet        -
     *    dls         -
     *    f2, f3      -
     *    pe          -
     *    pgh         -
     *    ph          -
     *    PInc        -
     *    pl          -
     *    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
     *    sll   , sls
     *    xls         -
     *    xnoh        -
     *    zf          -
     *    zm          -
     *
     *  coupling      :
     *    none.
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     * ----------------------------------------------------------------------------
     */
    static dpper_(options) {
        const { e3, ee2, peo, pgho, pho, PInco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos, } = options.satrec;
        let { ep, inclp, nodep, argpp, mp } = options;
        const { opsmode = _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.IMPROVED, init } = options;
        //  ---------------------- constants -----------------------------
        /** Ootk -- Some variables imported from outside the class at the top */
        const zns = 1.19459e-5;
        const zes = 0.01675;
        const znl = 1.5835218e-4;
        const zel = 0.0549;
        //  --------------- calculate time varying periodics -----------
        let zm = zmos + zns * t;
        // Be sure that the initial call has time set to zero
        if (init) {
            zm = zmos;
        }
        let zf = zm + 2.0 * zes * Math.sin(zm);
        let sinzf = Math.sin(zf);
        let f2 = 0.5 * sinzf * sinzf - 0.25;
        let f3 = -0.5 * sinzf * Math.cos(zf);
        const ses = se2 * f2 + se3 * f3;
        const sis = si2 * f2 + si3 * f3;
        const sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
        const sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
        const shs = sh2 * f2 + sh3 * f3;
        zm = zmol + znl * t;
        if (init) {
            zm = zmol;
        }
        zf = zm + 2.0 * zel * Math.sin(zm);
        sinzf = Math.sin(zf);
        f2 = 0.5 * sinzf * sinzf - 0.25;
        f3 = -0.5 * sinzf * Math.cos(zf);
        const sel = ee2 * f2 + e3 * f3;
        const sil = xi2 * f2 + xi3 * f3;
        const sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
        const sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
        const shll = xh2 * f2 + xh3 * f3;
        let pe = ses + sel;
        let PInc = sis + sil;
        let pl = sls + sll;
        let pgh = sghs + sghl;
        let ph = shs + shll;
        if (!init) {
            pe -= peo;
            PInc -= PInco;
            pl -= plo;
            pgh -= pgho;
            ph -= pho;
            inclp += PInc;
            ep += pe;
            const sinip = Math.sin(inclp);
            const cosip = Math.cos(inclp);
            /* ----------------- apply periodics directly ------------ */
            /*
             * Sgp4fix for lyddane choice
             * strn3 used original inclination - this is technically feasible
             * gsfc used perturbed inclination - also technically feasible
             * probably best to readjust the 0.2 limit value and limit discontinuity
             * 0.2 rad = 11.45916 deg
             * use next line for original strn3 approach and original inclination
             * if (inclo >= 0.2)
             * use next line for gsfc version and perturbed inclination
             */
            if (inclp >= 0.2) {
                ph /= sinip;
                pgh -= cosip * ph;
                argpp += pgh;
                nodep += ph;
                mp += pl;
            }
            else {
                //  ---- apply periodics with lyddane modification ----
                const sinop = Math.sin(nodep);
                const cosop = Math.cos(nodep);
                let alfdp = sinip * sinop;
                let betdp = sinip * cosop;
                const dalf = ph * cosop + PInc * cosip * sinop;
                const dbet = -ph * sinop + PInc * cosip * cosop;
                alfdp += dalf;
                betdp += dbet;
                nodep %= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                /*
                 *  Sgp4fix for afspc written intrinsic functions
                 *  nodep used without a trigonometric function ahead
                 */
                /* istanbul ignore next */
                if (nodep < 0.0 && opsmode === 'a') {
                    nodep += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                }
                let xls = mp + argpp + cosip * nodep;
                const dls = pl + pgh - PInc * nodep * sinip;
                xls += dls;
                const xnoh = nodep;
                nodep = Math.atan2(alfdp, betdp);
                /*
                 *  Sgp4fix for afspc written intrinsic functions
                 *  nodep used without a trigonometric function ahead
                 */
                /* istanbul ignore next */
                if (nodep < 0.0 && opsmode === 'a') {
                    nodep += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                }
                if (Math.abs(xnoh - nodep) > _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI) {
                    /* istanbul ignore next */
                    if (nodep < xnoh) {
                        nodep += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                    }
                    else {
                        nodep -= _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                    }
                }
                mp += pl;
                argpp = xls - mp - cosip * nodep;
            }
        } // If !init
        return {
            ep,
            inclp,
            nodep,
            argpp,
            mp,
        };
    }
    /*
     *-----------------------------------------------------------------------------
     *
     *                           procedure dscom
     *
     *  this procedure provides deep space common items used by both the secular
     *    and periodics subroutines.  input is provided as shown. this routine
     *    used to be called dpper, but the functions inside weren't well organized.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    epoch       -
     *    ep          - eccentricity
     *    argpp       - argument of perigee
     *    tc          -
     *    inclp       - inclination
     *    nodep       - right ascension of ascending node
     *    np          - mean motion
     *
     *  outputs       :
     *    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
     *    day         -
     *    e3          -
     *    ee2         -
     *    em          - eccentricity
     *    emsq        - eccentricity squared
     *    gam         -
     *    peo         -
     *    pgho        -
     *    pho         -
     *    PInco       -
     *    plo         -
     *    rtemsq      -
     *    se2, se3         -
     *    sgh2, sgh3, sgh4        -
     *    sh2, sh3, si2, si3, sl2, sl3, sl4         -
     *    s1, s2, s3, s4, s5, s6, s7          -
     *    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
     *    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
     *    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
     *    nm          - mean motion
     *    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
     *    zmol        -
     *    zmos        -
     *
     *  locals        :
     *    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
     *    betasq      -
     *    cc          -
     *    ctem, stem        -
     *    x1, x2, x3, x4, x5, x6, x7, x8          -
     *    xnodce      -
     *    xnoi        -
     *    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
     *    zcosi  , zsini  , zcosil , zsinil ,
     *    zx          -
     *    zy          -
     *
     *  coupling      :
     *    none.
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static dscom_(options) {
        const { epoch, ep, argpp, tc, inclp, nodep, np } = options;
        // -------------------------- constants -------------------------
        /** Ootk -- Some variables imported from outside the class at the top */
        const zes = 0.01675;
        const zel = 0.0549;
        const c1ss = 2.9864797e-6;
        const c1l = 4.7968065e-7;
        const zsinis = 0.39785416;
        const zcosis = 0.91744867;
        const zcosgs = 0.1945905;
        const zsings = -0.98088458;
        //  --------------------- local variables ------------------------
        let s1 = 0, s2 = 0, s3 = 0, s4 = 0, s5 = 0, s6 = 0, s7 = 0, ss1 = 0, ss2 = 0, ss3 = 0, ss4 = 0, ss5 = 0, ss6 = 0, ss7 = 0, sz1 = 0, sz11 = 0, sz12 = 0, sz13 = 0, sz2 = 0, sz21 = 0, sz22 = 0, sz23 = 0, sz3 = 0, sz31 = 0, sz32 = 0, sz33 = 0, z1 = 0, z11 = 0, z12 = 0, z13 = 0, z2 = 0, z21 = 0, z22 = 0, z23 = 0, z3 = 0, z31 = 0, z32 = 0, z33 = 0;
        const nm = np;
        const em = ep;
        const snodm = Math.sin(nodep);
        const cnodm = Math.cos(nodep);
        const sinomm = Math.sin(argpp);
        const cosomm = Math.cos(argpp);
        const sinim = Math.sin(inclp);
        const cosim = Math.cos(inclp);
        const emsq = em * em;
        const betasq = 1.0 - emsq;
        const rtemsq = Math.sqrt(betasq);
        //  ----------------- initialize lunar solar terms ---------------
        const peo = 0.0;
        const PInco = 0.0;
        const plo = 0.0;
        const pgho = 0.0;
        const pho = 0.0;
        const day = epoch + 18261.5 + tc / 1440.0;
        const xnodce = (4.523602 - 9.2422029e-4 * day) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        const stem = Math.sin(xnodce);
        const ctem = Math.cos(xnodce);
        const zcosil = 0.91375164 - 0.03568096 * ctem;
        const zsinil = Math.sqrt(1.0 - zcosil * zcosil);
        const zsinhl = (0.089683511 * stem) / zsinil;
        const zcoshl = Math.sqrt(1.0 - zsinhl * zsinhl);
        const gam = 5.8351514 + 0.001944368 * day;
        let zx = (0.39785416 * stem) / zsinil;
        const zy = zcoshl * ctem + 0.91744867 * zsinhl * stem;
        zx = Math.atan2(zx, zy);
        zx += gam - xnodce;
        const zcosgl = Math.cos(zx);
        const zsingl = Math.sin(zx);
        //  ------------------------- do solar terms ---------------------
        let zcosg = zcosgs;
        let zsing = zsings;
        let zcosi = zcosis;
        let zsini = zsinis;
        let zcosh = cnodm;
        let zsinh = snodm;
        let cc = c1ss;
        const xnoi = 1.0 / nm;
        for (let lsflg = 1; lsflg <= 2; lsflg++) {
            const a1 = zcosg * zcosh + zsing * zcosi * zsinh;
            const a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
            const a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
            const a8 = zsing * zsini;
            const a9 = zsing * zsinh + zcosg * zcosi * zcosh;
            const a10 = zcosg * zsini;
            const a2 = cosim * a7 + sinim * a8;
            const a4 = cosim * a9 + sinim * a10;
            const a5 = -sinim * a7 + cosim * a8;
            const a6 = -sinim * a9 + cosim * a10;
            const x1 = a1 * cosomm + a2 * sinomm;
            const x2 = a3 * cosomm + a4 * sinomm;
            const x3 = -a1 * sinomm + a2 * cosomm;
            const x4 = -a3 * sinomm + a4 * cosomm;
            const x5 = a5 * sinomm;
            const x6 = a6 * sinomm;
            const x7 = a5 * cosomm;
            const x8 = a6 * cosomm;
            z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
            z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
            z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
            z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
            z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
            z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
            z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
            z12 = -6.0 * (a1 * a6 + a3 * a5) + emsq * (-24.0 * (x2 * x7 + x1 * x8) + -6.0 * (x3 * x6 + x4 * x5));
            z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
            z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
            z22 = 6.0 * (a4 * a5 + a2 * a6) + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
            z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
            z1 = z1 + z1 + betasq * z31;
            z2 = z2 + z2 + betasq * z32;
            z3 = z3 + z3 + betasq * z33;
            s3 = cc * xnoi;
            s2 = (-0.5 * s3) / rtemsq;
            s4 = s3 * rtemsq;
            s1 = -15.0 * em * s4;
            s5 = x1 * x3 + x2 * x4;
            s6 = x2 * x3 + x1 * x4;
            s7 = x2 * x4 - x1 * x3;
            //  ----------------------- do lunar terms -------------------
            if (lsflg === 1) {
                ss1 = s1;
                ss2 = s2;
                ss3 = s3;
                ss4 = s4;
                ss5 = s5;
                ss6 = s6;
                ss7 = s7;
                sz1 = z1;
                sz2 = z2;
                sz3 = z3;
                sz11 = z11;
                sz12 = z12;
                sz13 = z13;
                sz21 = z21;
                sz22 = z22;
                sz23 = z23;
                sz31 = z31;
                sz32 = z32;
                sz33 = z33;
                zcosg = zcosgl;
                zsing = zsingl;
                zcosi = zcosil;
                zsini = zsinil;
                zcosh = zcoshl * cnodm + zsinhl * snodm;
                zsinh = snodm * zcoshl - cnodm * zsinhl;
                cc = c1l;
            }
        }
        const zmol = (4.7199672 + (0.2299715 * day - gam)) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        const zmos = (6.2565837 + 0.017201977 * day) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        //  ------------------------ do solar terms ----------------------
        const se2 = 2.0 * ss1 * ss6;
        const se3 = 2.0 * ss1 * ss7;
        const si2 = 2.0 * ss2 * sz12;
        const si3 = 2.0 * ss2 * (sz13 - sz11);
        const sl2 = -2.0 * ss3 * sz2;
        const sl3 = -2.0 * ss3 * (sz3 - sz1);
        const sl4 = -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
        const sgh2 = 2.0 * ss4 * sz32;
        const sgh3 = 2.0 * ss4 * (sz33 - sz31);
        const sgh4 = -18.0 * ss4 * zes;
        const sh2 = -2.0 * ss2 * sz22;
        const sh3 = -2.0 * ss2 * (sz23 - sz21);
        //  ------------------------ do lunar terms ----------------------
        const ee2 = 2.0 * s1 * s6;
        const e3 = 2.0 * s1 * s7;
        const xi2 = 2.0 * s2 * z12;
        const xi3 = 2.0 * s2 * (z13 - z11);
        const xl2 = -2.0 * s3 * z2;
        const xl3 = -2.0 * s3 * (z3 - z1);
        const xl4 = -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
        const xgh2 = 2.0 * s4 * z32;
        const xgh3 = 2.0 * s4 * (z33 - z31);
        const xgh4 = -18.0 * s4 * zel;
        const xh2 = -2.0 * s2 * z22;
        const xh3 = -2.0 * s2 * (z23 - z21);
        return {
            snodm,
            cnodm,
            sinim,
            cosim,
            sinomm,
            cosomm,
            day,
            e3,
            ee2,
            em,
            emsq,
            gam,
            peo,
            pgho,
            pho,
            PInco,
            plo,
            rtemsq,
            se2,
            se3,
            sgh2,
            sgh3,
            sgh4,
            sh2,
            sh3,
            si2,
            si3,
            sl2,
            sl3,
            sl4,
            s1,
            s2,
            s3,
            s4,
            s5,
            s6,
            s7,
            ss1,
            ss2,
            ss3,
            ss4,
            ss5,
            ss6,
            ss7,
            sz1,
            sz2,
            sz3,
            sz11,
            sz12,
            sz13,
            sz21,
            sz22,
            sz23,
            sz31,
            sz32,
            sz33,
            xgh2,
            xgh3,
            xgh4,
            xh2,
            xh3,
            xi2,
            xi3,
            xl2,
            xl3,
            xl4,
            nm,
            z1,
            z2,
            z3,
            z11,
            z12,
            z13,
            z21,
            z22,
            z23,
            z31,
            z32,
            z33,
            zmol,
            zmos,
        };
    }
    /*
     *-----------------------------------------------------------------------------
     *
     *                           procedure dsinit
     *
     *  this procedure provides deep space contributions to mean motion dot due
     *    to geopotential resonance with half day and one day orbits.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    cosim, sinim-
     *    emsq        - eccentricity squared
     *    argpo       - argument of perigee
     *    s1, s2, s3, s4, s5      -
     *    ss1, ss2, ss3, ss4, ss5 -
     *    sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
     *    t           - time
     *    tc          -
     *    gsto        - greenwich sidereal time                   rad
     *    mo          - mean anomaly
     *    mdot        - mean anomaly dot (rate)
     *    no          - mean motion
     *    nodeo       - right ascension of ascending node
     *    nodedot     - right ascension of ascending node dot (rate)
     *    xPIdot      -
     *    z1, z3, z11, z13, z21, z23, z31, z33 -
     *    eccm        - eccentricity
     *    argpm       - argument of perigee
     *    inclm       - inclination
     *    mm          - mean anomaly
     *    xn          - mean motion
     *    nodem       - right ascension of ascending node
     *
     *  outputs       :
     *    em          - eccentricity
     *    argpm       - argument of perigee
     *    inclm       - inclination
     *    mm          - mean anomaly
     *    nm          - mean motion
     *    nodem       - right ascension of ascending node
     *    irez        - flag for resonance           0-none, 1-one day, 2-half day
     *    atime       -
     *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433    -
     *    dedt        -
     *    didt        -
     *    dmdt        -
     *    dndt        -
     *    dnodt       -
     *    domdt       -
     *    del1, del2, del3        -
     *    ses  , sghl , sghs , sgs  , shl  , shs  , sis  , sls
     *    theta       -
     *    xfact       -
     *    xlamo       -
     *    xli         -
     *    xni
     *
     *  locals        :
     *    ainv2       -
     *    aonv        -
     *    cosisq      -
     *    eoc         -
     *    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543  -
     *    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533  -
     *    sini2       -
     *    temp        -
     *    temp1       -
     *    theta       -
     *    xno2        -
     *
     *  coupling      :
     *    getgravconst
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static dsinit_(options) {
        /*
         * Sgp4fix just send in xke as a constant and eliminate getgravconst call
         * gravconsttype whichconst,
         */
        const { xke, cosim, argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, t, tc, gsto, mo, mdot, no, nodeo, nodedot, xPIdot, z1, z3, z11, z13, z21, z23, z31, z33, ecco, eccsq, } = options;
        let { emsq, em, argpm, inclm, mm, nm, nodem, atime, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, del1, del2, del3, xfact, xlamo, xli, xni, } = options;
        /* --------------------- local variables ------------------------ */
        /** Ootk -- Some variables imported from outside the class at the top */
        const q22 = 1.7891679e-6;
        const q31 = 2.1460748e-6;
        const q33 = 2.2123015e-7;
        const root22 = 1.7891679e-6;
        const root44 = 7.3636953e-9;
        const root54 = 2.1765803e-9;
        const rptim = 4.37526908801129966e-3; // Equates to 7.29211514668855e-5 rad/sec
        const root32 = 3.7393792e-7;
        const root52 = 1.1428639e-7;
        const x2o3 = 2.0 / 3.0;
        const znl = 1.5835218e-4;
        const zns = 1.19459e-5;
        /*
         * Sgp4fix identify constants and allow alternate values
         * just xke is used here so pass it in rather than have multiple calls
         * getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
         */
        // -------------------- deep space initialization ------------
        let irez = 0;
        if (nm < 0.0052359877 && nm > 0.0034906585) {
            irez = 1;
        }
        if (nm >= 8.26e-3 && nm <= 9.24e-3 && em >= 0.5) {
            irez = 2;
        }
        // ------------------------ do solar terms -------------------
        const ses = ss1 * zns * ss5;
        const sis = ss2 * zns * (sz11 + sz13);
        const sls = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
        const sghs = ss4 * zns * (sz31 + sz33 - 6.0);
        let shs = -zns * ss2 * (sz21 + sz23);
        // Sgp4fix for 180 deg incl
        if (inclm < 5.2359877e-2 || inclm > _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI - 5.2359877e-2) {
            shs = 0.0;
        }
        if (sinim !== 0.0) {
            shs /= sinim;
        }
        const sgs = sghs - cosim * shs;
        // ------------------------- do lunar terms ------------------
        const dedt = ses + s1 * znl * s5;
        const didt = sis + s2 * znl * (z11 + z13);
        const dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
        const sghl = s4 * znl * (z31 + z33 - 6.0);
        let shll = -znl * s2 * (z21 + z23);
        // Sgp4fix for 180 deg incl
        if (inclm < 5.2359877e-2 || inclm > _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.PI - 5.2359877e-2) {
            shll = 0.0;
        }
        let domdt = sgs + sghl;
        let dnodt = shs;
        if (sinim !== 0.0) {
            domdt -= (cosim / sinim) * shll;
            dnodt += shll / sinim;
        }
        // ----------- calculate deep space resonance effects --------
        const dndt = 0.0;
        const theta = (gsto + tc * rptim) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        em += dedt * t;
        inclm += didt * t;
        argpm += domdt * t;
        nodem += dnodt * t;
        mm += dmdt * t;
        /*
         * Sgp4fix for negative inclinations
         * the following if statement should be commented out
         * if (inclm < 0.0)
         * {
         *   inclm  = -inclm;
         *   argpm  = argpm - PI;
         *   nodem = nodem + PI;
         * }
         */
        // -------------- initialize the resonance terms -------------
        if (irez !== 0) {
            const aonv = (nm / xke) ** x2o3;
            // ---------- geopotential resonance for 12 hour orbits ------
            if (irez === 2) {
                const cosisq = cosim * cosim;
                const emo = em;
                em = ecco;
                const emsqo = emsq;
                emsq = eccsq;
                const eoc = em * emsq;
                const g201 = -0.306 - (em - 0.64) * 0.44;
                let g211, g310, g322, g410, g422, g520, g521, g532, g533;
                if (em <= 0.65) {
                    g211 = 3.616 - 13.247 * em + 16.29 * emsq;
                    g310 = -19.302 + 117.39 * em - 228.419 * emsq + 156.591 * eoc;
                    g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
                    g410 = -41.122 + 242.694 * em - 471.094 * emsq + 313.953 * eoc;
                    g422 = -146.407 + 841.88 * em - 1629.014 * emsq + 1083.435 * eoc;
                    g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.276 * eoc;
                }
                else {
                    g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
                    g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
                    g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
                    g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
                    g422 = -3581.69 + 16178.11 * em - 24462.77 * emsq + 12422.52 * eoc;
                    if (em > 0.715) {
                        g520 = -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
                    }
                    else {
                        g520 = 1464.74 - 4664.75 * em + 3763.64 * emsq;
                    }
                }
                if (em < 0.7) {
                    g533 = -919.2277 + 4988.61 * em - 9064.77 * emsq + 5542.21 * eoc;
                    g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                    g532 = -853.666 + 4690.25 * em - 8624.77 * emsq + 5341.4 * eoc;
                }
                else {
                    g533 = -37995.78 + 161616.52 * em - 229838.2 * emsq + 109377.94 * eoc;
                    g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                    g532 = -40023.88 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
                }
                const sini2 = sinim * sinim;
                const f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
                const f221 = 1.5 * sini2;
                const f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
                const f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
                const f441 = 35.0 * sini2 * f220;
                const f442 = 39.375 * sini2 * sini2;
                const f522 = 9.84375 *
                    sinim *
                    (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq) + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
                const f523 = sinim *
                    (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
                const f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
                const f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
                const xno2 = nm * nm;
                const ainv2 = aonv * aonv;
                let temp1 = 3.0 * xno2 * ainv2;
                let temp = temp1 * root22;
                d2201 = temp * f220 * g201;
                d2211 = temp * f221 * g211;
                temp1 *= aonv;
                temp = temp1 * root32;
                d3210 = temp * f321 * g310;
                d3222 = temp * f322 * g322;
                temp1 *= aonv;
                temp = 2.0 * temp1 * root44;
                d4410 = temp * f441 * g410;
                d4422 = temp * f442 * g422;
                temp1 *= aonv;
                temp = temp1 * root52;
                d5220 = temp * f522 * g520;
                d5232 = temp * f523 * g532;
                temp = 2.0 * temp1 * root54;
                d5421 = temp * f542 * g521;
                d5433 = temp * f543 * g533;
                xlamo = (mo + nodeo + nodeo - (theta + theta)) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                xfact = mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
                em = emo;
                emsq = emsqo;
            }
            //  ---------------- synchronous resonance terms --------------
            if (irez === 1) {
                const g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
                const g310 = 1.0 + 2.0 * emsq;
                const g300 = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
                const f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                const f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
                let f330 = 1.0 + cosim;
                f330 *= 1.875 * f330 * f330;
                del1 = 3.0 * nm * nm * aonv * aonv;
                del2 = 2.0 * del1 * f220 * g200 * q22;
                del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
                del1 = del1 * f311 * g310 * q31 * aonv;
                xlamo = (mo + nodeo + argpo - theta) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
                xfact = mdot + xPIdot + dmdt + domdt + dnodt - (no + rptim);
            }
            //  ------------ for sgp4, initialize the integrator ----------
            xli = xlamo;
            xni = no;
            atime = 0.0;
            nm = no + dndt;
        }
        return {
            em,
            argpm,
            inclm,
            mm,
            nm,
            nodem,
            irez,
            atime,
            d2201,
            d2211,
            d3210,
            d3222,
            d4410,
            d4422,
            d5220,
            d5232,
            d5421,
            d5433,
            dedt,
            didt,
            dmdt,
            dndt,
            dnodt,
            domdt,
            del1,
            del2,
            del3,
            xfact,
            xlamo,
            xli,
            xni,
        };
    }
    /*
     *-----------------------------------------------------------------------------
     *
     *                           procedure dspace
     *
     *  this procedure provides deep space contributions to mean elements for
     *    perturbing third body.  these effects have been averaged over one
     *    revolution of the sun and moon.  for earth resonance effects, the
     *    effects have been averaged over no revolutions of the satellite.
     *    (mean motion)
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
     *    dedt        -
     *    del1, del2, del3  -
     *    didt        -
     *    dmdt        -
     *    dnodt       -
     *    domdt       -
     *    irez        - flag for resonance           0-none, 1-one day, 2-half day
     *    argpo       - argument of perigee
     *    argpdot     - argument of perigee dot (rate)
     *    t           - time
     *    tc          -
     *    gsto        - gst
     *    xfact       -
     *    xlamo       -
     *    no          - mean motion
     *    atime       -
     *    em          - eccentricity
     *    ft          -
     *    argpm       - argument of perigee
     *    inclm       - inclination
     *    xli         -
     *    mm          - mean anomaly
     *    xni         - mean motion
     *    nodem       - right ascension of ascending node
     *
     *  outputs       :
     *    atime       -
     *    em          - eccentricity
     *    argpm       - argument of perigee
     *    inclm       - inclination
     *    xli         -
     *    mm          - mean anomaly
     *    xni         -
     *    nodem       - right ascension of ascending node
     *    dndt        -
     *    nm          - mean motion
     *
     *  locals        :
     *    delt        -
     *    ft          -
     *    theta       -
     *    x2li        -
     *    x2omi       -
     *    xl          -
     *    xldot       -
     *    xnddt       -
     *    xndt        -
     *    xomi        -
     *
     *  coupling      :
     *    none        -
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static dspace_(em, argpm, inclm, mm, nodem, nm, satrec) {
        let { atime, xli, xni, } = satrec;
        const { dedt, didt, dmdt, dnodt, domdt, irez, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, xfact, xlamo, gsto, argpo, t, no, argpdot, del1, del2, del3, } = satrec;
        //  ----------- calculate deep space resonance effects -----------
        let dndt = 0.0;
        const theta = (gsto + t * rptim) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
        // Apply time-dependent perturbations
        em += dedt * t;
        inclm += didt * t;
        argpm += domdt * t;
        nodem += dnodt * t;
        mm += dmdt * t;
        if (irez !== 0) {
            // Simplified deep space resonance handling
            if (atime === 0.0 || t * atime <= 0.0 || Math.abs(t) < Math.abs(atime)) {
                atime = 0.0;
                xni = no;
                xli = xlamo;
            }
            const delt = t > 0.0 ? stepp : stepn;
            let ft = 0;
            let x2li = 0;
            let x2omi = 0;
            let xldot = 0;
            let xnddt = 0;
            let xndt = 0;
            let xomi = 0;
            let iretn = true; // Added for do loop
            while (iretn) {
                /*
                 *  ------------------- dot terms calculated -------------
                 *  ----------- near - synchronous resonance terms -------
                 */
                if (irez !== 2) {
                    xndt =
                        del1 * Math.sin(xli - fasx2) + del2 * Math.sin(2.0 * (xli - fasx4)) + del3 * Math.sin(3.0 * (xli - fasx6));
                    xldot = xni + xfact;
                    xnddt =
                        del1 * Math.cos(xli - fasx2) +
                            2.0 * del2 * Math.cos(2.0 * (xli - fasx4)) +
                            3.0 * del3 * Math.cos(3.0 * (xli - fasx6));
                    xnddt *= xldot;
                }
                else {
                    // --------- near - half-day resonance terms --------
                    xomi = argpo + argpdot * atime;
                    x2omi = xomi + xomi;
                    x2li = xli + xli;
                    xndt =
                        d2201 * Math.sin(x2omi + xli - g22) +
                            d2211 * Math.sin(xli - g22) +
                            d3210 * Math.sin(xomi + xli - g32) +
                            d3222 * Math.sin(-xomi + xli - g32) +
                            d4410 * Math.sin(x2omi + x2li - g44) +
                            d4422 * Math.sin(x2li - g44) +
                            d5220 * Math.sin(xomi + xli - g52) +
                            d5232 * Math.sin(-xomi + xli - g52) +
                            d5421 * Math.sin(xomi + x2li - g54) +
                            d5433 * Math.sin(-xomi + x2li - g54);
                    xldot = xni + xfact;
                    xnddt =
                        d2201 * Math.cos(x2omi + xli - g22) +
                            d2211 * Math.cos(xli - g22) +
                            d3210 * Math.cos(xomi + xli - g32) +
                            d3222 * Math.cos(-xomi + xli - g32) +
                            d5220 * Math.cos(xomi + xli - g52) +
                            d5232 * Math.cos(-xomi + xli - g52) +
                            2.0 *
                                (d4410 * Math.cos(x2omi + x2li - g44) +
                                    d4422 * Math.cos(x2li - g44) +
                                    d5421 * Math.cos(xomi + x2li - g54) +
                                    d5433 * Math.cos(-xomi + x2li - g54));
                    xnddt *= xldot;
                }
                /*
                 *  ----------------------- integrator -------------------
                 *  sgp4fix move end checks to end of routine
                 */
                if (Math.abs(t - atime) < stepp) {
                    ft = t - atime;
                    iretn = false;
                }
                else {
                    xli += xldot * delt + xndt * step2;
                    xni += xndt * delt + xnddt * step2;
                    atime += delt;
                }
            }
            nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
            const xl = xli + xldot * ft + xndt * ft * ft * 0.5;
            if (irez !== 1) {
                mm = xl - 2.0 * nodem + 2.0 * theta;
                dndt = nm - no;
            }
            else {
                mm = xl - nodem - argpm + theta;
                dndt = nm - no;
            }
            nm = no + dndt;
        }
        return [em, argpm, inclm, mm, nodem, nm];
    }
    /*
     * -----------------------------------------------------------------------------
     *
     *                           function getgravconst
     *
     *  this function gets constants for the propagator. note that mu is identified to
     *    facilitiate comparisons with newer models. the common useage is wgs72.
     *
     *  author        : david vallado                  719-573-2600   21 jul 2006
     *
     *  inputs        :
     *    whichconst  - which set of constants to use  wgs72old, wgs72, wgs84
     *
     *  outputs       :
     *    tumin       - minutes in one time unit
     *    mu          - earth gravitational parameter
     *    radiusearthkm - radius of the earth in km
     *    xke         - reciprocal of tumin
     *    j2, j3, j4  - un-normalized zonal harmonic values
     *    j3oj2       - j3 divided by j2
     *
     *  locals        :
     *
     *  coupling      :
     *    none
     *
     *  references    :
     *    norad spacetrack report #3
     *    vallado, crawford, hujsak, kelso  2006
     * ---------------------------------------------------------------------------
     */
    static getgravconst_(whichconst) {
        let j2, j3, j3oj2, j4, mus, radiusearthkm, tumin, xke;
        switch (whichconst) {
            // -- wgs-72 low precision str#3 constants --
            case 'wgs72old':
                mus = 398600.79964; // In km3 / s2
                radiusearthkm = 6378.135; // Km
                xke = 0.0743669161; // Reciprocal of tumin
                tumin = 1.0 / xke;
                j2 = 0.001082616;
                j3 = -0.00000253881;
                j4 = -0.00000165597;
                j3oj2 = j3 / j2;
                break;
            // ------------ wgs-72 constants ------------
            case 'wgs72':
                mus = 398600.8; // In km3 / s2
                radiusearthkm = 6378.135; // Km
                xke = 60.0 / Math.sqrt((radiusearthkm * radiusearthkm * radiusearthkm) / mus);
                tumin = 1.0 / xke;
                j2 = 0.001082616;
                j3 = -0.00000253881;
                j4 = -0.00000165597;
                j3oj2 = j3 / j2;
                break;
            case 'wgs84':
                // ------------ wgs-84 constants ------------
                mus = 398600.5; // In km3 / s2
                radiusearthkm = 6378.137; // Km
                xke = 60.0 / Math.sqrt((radiusearthkm * radiusearthkm * radiusearthkm) / mus);
                tumin = 1.0 / xke;
                j2 = 0.00108262998905;
                j3 = -0.00000253215306;
                j4 = -0.00000161098761;
                j3oj2 = j3 / j2;
                break;
            default:
                throw new Error(`unknown gravity option ${whichconst}`);
        }
        return {
            tumin,
            mus,
            radiusearthkm,
            xke,
            j2,
            j3,
            j4,
            j3oj2,
        };
    }
    /*
     *-----------------------------------------------------------------------------
     *
     *                           procedure initl
     *
     *  this procedure initializes the sgp4 propagator. all the initialization is
     *    consolidated here instead of having multiple loops inside other routines.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    satn        - satellite number - not needed, placed in satrec
     *    xke         - reciprocal of tumin
     *    j2          - j2 zonal harmonic
     *    ecco        - eccentricity                           0.0 - 1.0
     *    epoch       - epoch time in days from jan 0, 1950. 0 hr
     *    inclo       - inclination of satellite
     *    no          - mean motion of satellite
     *
     *  outputs       :
     *    ainv        - 1.0 / a
     *    ao          - semi major axis
     *    con41       -
     *    con42       - 1.0 - 5.0 cos(i)
     *    cosio       - cosine of inclination
     *    cosio2      - cosio squared
     *    eccsq       - eccentricity squared
     *    method      - flag for deep space                    'd', 'n'
     *    omeosq      - 1.0 - ecco * ecco
     *    posq        - semi-parameter squared
     *    rp          - radius of perigee
     *    rteosq      - square root of (1.0 - ecco*ecco)
     *    sinio       - sine of inclination
     *    gsto        - gst at time of observation               rad
     *    no          - mean motion of satellite
     *
     *  locals        :
     *    ak          -
     *    d1          -
     *    del         -
     *    adel        -
     *    po          -
     *
     *  coupling      :
     *    getgravconst- no longer used
     *    gstime      - find greenwich sidereal time from the julian date
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static initl_(satrec, epoch) {
        /*
         * Sgp4fix satn not needed. include in satrec in case needed later
         * int satn,
         * sgp4fix just pass in xke and j2
         * gravconsttype whichconst,
         */
        const { operationmode, ecco, inclo, xke, j2 } = satrec;
        let { no } = satrec;
        /* --------------------- local variables ------------------------ */
        const x2o3 = 2.0 / 3.0;
        // Sgp4fix use old way of finding gst
        /** Ootk -- Some variables imported from outside the class at the top */
        /*
         * ----------------------- earth constants ---------------------
         * sgp4fix identify constants and allow alternate values
         * only xke and j2 are used here so pass them in directly
         * getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 )
         */
        // ------------- calculate auxillary epoch quantities ----------
        const eccsq = ecco * ecco;
        const omeosq = 1.0 - eccsq;
        const rteosq = Math.sqrt(omeosq);
        const cosio = Math.cos(inclo);
        const cosio2 = cosio * cosio;
        // ------------------ un-kozai the mean motion -----------------
        const ak = (xke / no) ** x2o3;
        const d1 = (0.75 * j2 * (3.0 * cosio2 - 1.0)) / (rteosq * omeosq);
        let delPrime = d1 / (ak * ak);
        const adel = ak * (1.0 - delPrime * delPrime - delPrime * (1.0 / 3.0 + (134.0 * delPrime * delPrime) / 81.0));
        delPrime = d1 / (adel * adel);
        no /= 1.0 + delPrime;
        const ao = (xke / no) ** x2o3;
        const sinio = Math.sin(inclo);
        const po = ao * omeosq;
        const con42 = 1.0 - 5.0 * cosio2;
        const con41 = -con42 - cosio2 - cosio2;
        const ainv = 1.0 / ao;
        const posq = po * po;
        const rp = ao * (1.0 - ecco);
        //  Sgp4fix modern approach to finding sidereal time
        /** Ootk -- Continue allowing AFSPC mode for SGP4 Validation */
        let gsto;
        if (operationmode === _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.AFSPC) {
            /*
             *  Sgp4fix use old way of finding gst
             *  count integer number of days from 0 jan 1970
             */
            const ts70 = epoch - 7305.0;
            const ds70 = Math.floor(ts70 + 1.0e-8);
            const tfrac = ts70 - ds70;
            //  Find greenwich location at epoch
            const c1 = 1.72027916940703639e-2;
            const thgr70 = 1.7321343856509374;
            const fk5r = 5.07551419432269442e-15;
            const c1p2p = c1 + _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
            gsto = (thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
            /* istanbul ignore next | This is no longer possible*/
            if (gsto < 0.0) {
                gsto += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
            }
        }
        else {
            const jdut1 = epoch + 2433281.5;
            const tut1 = (jdut1 - 2451545.0) / 36525.0;
            gsto = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
                (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841; // Sec
            gsto = ((gsto * _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.DEG2RAD) / 240.0) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU; // 360/86400 = 1/240, to deg, to rad
            //  ------------------------ check quadrants ---------------------
            if (gsto < 0.0) {
                gsto += _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU;
            }
        }
        return {
            no,
            ainv,
            ao,
            con41,
            con42,
            cosio,
            cosio2,
            eccsq,
            omeosq,
            posq,
            rp,
            rteosq,
            sinio,
            gsto,
        };
    }
    /*
     *-----------------------------------------------------------------------------
     *
     *                             procedure sgp4init
     *
     *  this procedure initializes variables for sgp4.
     *
     *  author        : david vallado                  719-573-2600   28 jun 2005
     *
     *  inputs        :
     *    opsmode     - mode of operation afspc or improved 'a', 'i'
     *    whichconst  - which set of constants to use  72, 84
     *    satn        - satellite number
     *    bstar       - sgp4 type drag coefficient              kg/m2er
     *    ecco        - eccentricity
     *    epoch       - epoch time in days from jan 0, 1950. 0 hr
     *    argpo       - argument of perigee (output if ds)
     *    inclo       - inclination
     *    mo          - mean anomaly (output if ds)
     *    no          - mean motion
     *    nodeo       - right ascension of ascending node
     *
     *  outputs       :
     *    rec      - common values for subsequent calls
     *    return code - non-zero on error.
     *                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
     *                   2 - mean motion less than 0.0
     *                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
     *                   4 - semi-latus rectum < 0.0
     *                   5 - epoch elements are sub-orbital
     *                   6 - satellite has decayed
     *
     *  locals        :
     *    cnodm  , snodm  , cosim  , sinim  , cosomm , sinomm
     *    cc1sq  , cc2    , cc3
     *    coef   , coef1
     *    cosio4      -
     *    day         -
     *    dndt        -
     *    em          - eccentricity
     *    emsq        - eccentricity squared
     *    eeta        -
     *    etasq       -
     *    gam         -
     *    argpm       - argument of perigee
     *    nodem       -
     *    inclm       - inclination
     *    mm          - mean anomaly
     *    nm          - mean motion
     *    perige      - perigee
     *    PInvsq      -
     *    psisq       -
     *    qzms24      -
     *    rtemsq      -
     *    s1, s2, s3, s4, s5, s6, s7          -
     *    sfour       -
     *    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
     *    sz1, sz2, sz3
     *    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
     *    tc          -
     *    temp        -
     *    temp1, temp2, temp3       -
     *    tsi         -
     *    xPIdot      -
     *    xhdot1      -
     *    z1, z2, z3          -
     *    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
     *
     *  coupling      :
     *    getgravconst-
     *    initl       -
     *    dscom       -
     *    dpper       -
     *    dsinit      -
     *    sgp4        -
     *
     *  references    :
     *    hoots, roehrich, norad spacetrack report #3 1980
     *    hoots, norad spacetrack report #6 1986
     *    hoots, schumacher and glover 2004
     *    vallado, crawford, hujsak, kelso  2006
     *----------------------------------------------------------------------------
     */
    static sgp4init_(satrec, options) {
        const { whichconst = Sgp4GravConstants.wgs72, opsmode = _enums_Sgp4OpsMode_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4OpsMode.IMPROVED, satn = satrec.satnum, epoch, xbstar, xecco, xargpo, xinclo, xndot, xnddot, xmo, xno, xnodeo, } = options;
        /* ------------------------ initialization --------------------- */
        /*
         * Sgp4fix divisor for divide by zero check on inclination
         * the old check used 1.0 + Math.cos(PI-1.0e-9), but then compared it to
         * 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
         */
        // ----------- set all near earth variables to zero ------------
        satrec.isimp = false;
        satrec.method = Sgp4Methods.NEAR_EARTH;
        satrec.aycof = 0.0;
        satrec.con41 = 0.0;
        satrec.cc1 = 0.0;
        satrec.cc4 = 0.0;
        satrec.cc5 = 0.0;
        satrec.d2 = 0.0;
        satrec.d3 = 0.0;
        satrec.d4 = 0.0;
        satrec.delmo = 0.0;
        satrec.eta = 0.0;
        satrec.argpdot = 0.0;
        satrec.omgcof = 0.0;
        satrec.sinmao = 0.0;
        satrec.t = 0.0;
        satrec.t2cof = 0.0;
        satrec.t3cof = 0.0;
        satrec.t4cof = 0.0;
        satrec.t5cof = 0.0;
        satrec.x1mth2 = 0.0;
        satrec.x7thm1 = 0.0;
        satrec.mdot = 0.0;
        satrec.nodedot = 0.0;
        satrec.xlcof = 0.0;
        satrec.xmcof = 0.0;
        satrec.nodecf = 0.0;
        // ----------- set all deep space variables to zero ------------
        satrec.irez = 0;
        satrec.d2201 = 0.0;
        satrec.d2211 = 0.0;
        satrec.d3210 = 0.0;
        satrec.d3222 = 0.0;
        satrec.d4410 = 0.0;
        satrec.d4422 = 0.0;
        satrec.d5220 = 0.0;
        satrec.d5232 = 0.0;
        satrec.d5421 = 0.0;
        satrec.d5433 = 0.0;
        satrec.dedt = 0.0;
        satrec.del1 = 0.0;
        satrec.del2 = 0.0;
        satrec.del3 = 0.0;
        satrec.didt = 0.0;
        satrec.dmdt = 0.0;
        satrec.dnodt = 0.0;
        satrec.domdt = 0.0;
        satrec.e3 = 0.0;
        satrec.ee2 = 0.0;
        satrec.peo = 0.0;
        satrec.pgho = 0.0;
        satrec.pho = 0.0;
        satrec.PInco = 0.0;
        satrec.plo = 0.0;
        satrec.se2 = 0.0;
        satrec.se3 = 0.0;
        satrec.sgh2 = 0.0;
        satrec.sgh3 = 0.0;
        satrec.sgh4 = 0.0;
        satrec.sh2 = 0.0;
        satrec.sh3 = 0.0;
        satrec.si2 = 0.0;
        satrec.si3 = 0.0;
        satrec.sl2 = 0.0;
        satrec.sl3 = 0.0;
        satrec.sl4 = 0.0;
        satrec.gsto = 0.0;
        satrec.xfact = 0.0;
        satrec.xgh2 = 0.0;
        satrec.xgh3 = 0.0;
        satrec.xgh4 = 0.0;
        satrec.xh2 = 0.0;
        satrec.xh3 = 0.0;
        satrec.xi2 = 0.0;
        satrec.xi3 = 0.0;
        satrec.xl2 = 0.0;
        satrec.xl3 = 0.0;
        satrec.xl4 = 0.0;
        satrec.xlamo = 0.0;
        satrec.zmol = 0.0;
        satrec.zmos = 0.0;
        satrec.atime = 0.0;
        satrec.xli = 0.0;
        satrec.xni = 0.0;
        /* ------------------------ earth constants ----------------------- */
        /*
         * Sgp4fix identify constants and allow alternate values
         * this is now the only call for the constants
         */
        const gravResults = Sgp4.getgravconst_(whichconst);
        satrec.tumin = gravResults.tumin;
        satrec.mus = gravResults.mus;
        satrec.radiusearthkm = gravResults.radiusearthkm;
        satrec.xke = gravResults.xke;
        satrec.j2 = gravResults.j2;
        satrec.j3 = gravResults.j3;
        satrec.j4 = gravResults.j4;
        satrec.j3oj2 = gravResults.j3oj2;
        satrec.vkmpersec = (satrec.radiusearthkm * satrec.xke) / 60.0;
        const { j2 } = gravResults;
        const { j4 } = gravResults;
        const { xke } = gravResults;
        const { j3oj2 } = gravResults;
        // -------------------------------------------------------------------------
        satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR;
        satrec.operationmode = opsmode;
        // New alpha5 or 9-digit number
        /**
         * Ootk -- Using JS code for string manipulation but same effect
         * Ex. A2525 = 102525
         * Ex. Z1234 = 351234
         */
        satrec.satnum = _main_js__WEBPACK_IMPORTED_MODULE_1__.Tle.convertA5to6Digit(satn);
        /*
         * Sgp4fix - note the following variables are also passed directly via satrec.
         * it is possible to streamline the sgp4init call by deleting the "x"
         * variables, but the user would need to set the satrec.* values first. we
         * include the additional assignments in case twoline2rv is not used.
         */
        satrec.bstar = xbstar;
        // Sgp4fix allow additional parameters in the struct
        satrec.ndot = xndot;
        satrec.nddot = xnddot;
        satrec.ecco = xecco;
        satrec.argpo = xargpo;
        satrec.inclo = xinclo;
        satrec.mo = xmo;
        // Sgp4fix rename variables to clarify which mean motion is intended
        satrec.no = xno;
        satrec.nodeo = xnodeo;
        /*
         * Single averaged mean elements
         * satrec.am = satrec.em = satrec.im = satrec.Om = satrec.mm = satrec.nm = 0.0;
         */
        /*
         * ------------------------ earth constants -----------------------
         * sgp4fix identify constants and allow alternate values
         * getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 );
         */
        const ss = 78.0 / satrec.radiusearthkm + 1.0;
        // Sgp4fix use multiply for speed instead of pow
        const qzms2ttemp = (120.0 - 78.0) / satrec.radiusearthkm;
        const qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
        satrec.init = true;
        satrec.t = 0.0;
        // Sgp4fix remove satn as it is not needed in initl
        const initlResult = Sgp4.initl_(satrec, epoch);
        const { ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio } = initlResult;
        satrec.no = initlResult.no;
        satrec.con41 = initlResult.con41;
        satrec.gsto = initlResult.gsto;
        satrec.a = (satrec.no * satrec.tumin) ** (-2.0 / 3.0);
        satrec.alta = satrec.a * (1.0 + satrec.ecco) - 1.0;
        satrec.altp = satrec.a * (1.0 - satrec.ecco) - 1.0;
        satrec.error = _sgp4_error_js__WEBPACK_IMPORTED_MODULE_3__.Sgp4ErrorCode.NO_ERROR;
        /*
         * Sgp4fix remove this check as it is unnecessary
         * the mrt check in sgp4 handles decaying satellite cases even if the starting
         * condition is below the surface of te earth
         * if (rp < 1.0)
         * {
         *   printf("// *** satn%d epoch elts sub-orbital ***\n", satn);
         *   satrec.error = 5;
         * }
         */
        if (omeosq >= 0.0 || satrec.no >= 0.0) {
            satrec.isimp = false;
            if (rp < 220.0 / satrec.radiusearthkm + 1.0) {
                satrec.isimp = true;
            }
            let sfour = ss;
            let qzms24 = qzms2t;
            const perigee = (rp - 1.0) * satrec.radiusearthkm;
            // - for perigees below 156 km, s and qoms2t are altered -
            if (perigee < 156.0) {
                sfour = perigee - 78.0;
                if (perigee < 98.0) {
                    sfour = 20.0;
                }
                // Sgp4fix use multiply for speed instead of pow
                const qzms24temp = (120.0 - sfour) / satrec.radiusearthkm;
                qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                sfour = sfour / satrec.radiusearthkm + 1.0;
            }
            const PInvsq = 1.0 / posq;
            const tsi = 1.0 / (ao - sfour);
            satrec.eta = ao * satrec.ecco * tsi;
            const etasq = satrec.eta * satrec.eta;
            const eeta = satrec.ecco * satrec.eta;
            const psisq = Math.abs(1.0 - etasq);
            const coef = qzms24 * tsi ** 4.0;
            const coef1 = coef / psisq ** 3.5;
            const cc2 = coef1 *
                satrec.no *
                (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
                    ((0.375 * j2 * tsi) / psisq) * satrec.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
            satrec.cc1 = satrec.bstar * cc2;
            let cc3 = 0.0;
            if (satrec.ecco > 1.0e-4) {
                cc3 = (-2.0 * coef * tsi * j3oj2 * satrec.no * sinio) / satrec.ecco;
            }
            satrec.x1mth2 = 1.0 - cosio2;
            satrec.cc4 =
                2.0 *
                    satrec.no *
                    coef1 *
                    ao *
                    omeosq *
                    (satrec.eta * (2.0 + 0.5 * etasq) +
                        satrec.ecco * (0.5 + 2.0 * etasq) -
                        ((j2 * tsi) / (ao * psisq)) *
                            (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) +
                                0.75 * satrec.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * Math.cos(2.0 * satrec.argpo)));
            satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
            const cosio4 = cosio2 * cosio2;
            const temp1 = 1.5 * j2 * PInvsq * satrec.no;
            const temp2 = 0.5 * temp1 * j2 * PInvsq;
            const temp3 = -0.46875 * j4 * PInvsq * PInvsq * satrec.no;
            satrec.mdot =
                satrec.no +
                    0.5 * temp1 * rteosq * satrec.con41 +
                    0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
            satrec.argpdot =
                -0.5 * temp1 * con42 +
                    0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                    temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
            const xhdot1 = -temp1 * cosio;
            satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
            const xPIdot = satrec.argpdot + satrec.nodedot;
            satrec.omgcof = satrec.bstar * cc3 * Math.cos(satrec.argpo);
            satrec.xmcof = 0.0;
            if (satrec.ecco > 1.0e-4) {
                satrec.xmcof = (-_utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.x2o3 * coef * satrec.bstar) / eeta;
            }
            satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
            satrec.t2cof = 1.5 * satrec.cc1;
            // Sgp4fix for divide by zero with xinco = 180 deg
            if (Math.abs(cosio + 1.0) > 1.5e-12) {
                satrec.xlcof = (-0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio)) / (1.0 + cosio);
            }
            else {
                satrec.xlcof = (-0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio)) / _utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.temp4;
            }
            satrec.aycof = -0.5 * j3oj2 * sinio;
            // Sgp4fix use multiply for speed instead of pow
            const delmotemp = 1.0 + satrec.eta * Math.cos(satrec.mo);
            satrec.delmo = delmotemp * delmotemp * delmotemp;
            satrec.sinmao = Math.sin(satrec.mo);
            satrec.x7thm1 = 7.0 * cosio2 - 1.0;
            // --------------- deep space initialization -------------
            if (_utils_constants_js__WEBPACK_IMPORTED_MODULE_2__.TAU / satrec.no >= 225.0) {
                satrec.method = Sgp4Methods.DEEP_SPACE;
                satrec.isimp = true;
                const tc = 0.0;
                const inclm = satrec.inclo;
                const dscomOptions = {
                    epoch,
                    ep: satrec.ecco,
                    argpp: satrec.argpo,
                    tc,
                    inclp: satrec.inclo,
                    nodep: satrec.nodeo,
                    np: satrec.no,
                    e3: satrec.e3,
                    ee2: satrec.ee2,
                    peo: satrec.peo,
                    pgho: satrec.pgho,
                    pho: satrec.pho,
                    PInco: satrec.PInco,
                    plo: satrec.plo,
                    se2: satrec.se2,
                    se3: satrec.se3,
                    sgh2: satrec.sgh2,
                    sgh3: satrec.sgh3,
                    sgh4: satrec.sgh4,
                    sh2: satrec.sh2,
                    sh3: satrec.sh3,
                    si2: satrec.si2,
                    si3: satrec.si3,
                    sl2: satrec.sl2,
                    sl3: satrec.sl3,
                    sl4: satrec.sl4,
                    xgh2: satrec.xgh2,
                    xgh3: satrec.xgh3,
                    xgh4: satrec.xgh4,
                    xh2: satrec.xh2,
                    xh3: satrec.xh3,
                    xi2: satrec.xi2,
                    xi3: satrec.xi3,
                    xl2: satrec.xl2,
                    xl3: satrec.xl3,
                    xl4: satrec.xl4,
                    zmol: satrec.zmol,
                    zmos: satrec.zmos,
                };
                const dscomResult = Sgp4.dscom_(dscomOptions);
                satrec.e3 = dscomResult.e3;
                satrec.ee2 = dscomResult.ee2;
                satrec.peo = dscomResult.peo;
                satrec.pgho = dscomResult.pgho;
                satrec.pho = dscomResult.pho;
                satrec.PInco = dscomResult.PInco;
                satrec.plo = dscomResult.plo;
                satrec.se2 = dscomResult.se2;
                satrec.se3 = dscomResult.se3;
                satrec.sgh2 = dscomResult.sgh2;
                satrec.sgh3 = dscomResult.sgh3;
                satrec.sgh4 = dscomResult.sgh4;
                satrec.sh2 = dscomResult.sh2;
                satrec.sh3 = dscomResult.sh3;
                satrec.si2 = dscomResult.si2;
                satrec.si3 = dscomResult.si3;
                satrec.sl2 = dscomResult.sl2;
                satrec.sl3 = dscomResult.sl3;
                satrec.sl4 = dscomResult.sl4;
                const { sinim, cosim, em, emsq, s1, s2, s3, s4, s5, ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, } = dscomResult;
                satrec.xgh2 = dscomResult.xgh2;
                satrec.xgh3 = dscomResult.xgh3;
                satrec.xgh4 = dscomResult.xgh4;
                satrec.xh2 = dscomResult.xh2;
                satrec.xh3 = dscomResult.xh3;
                satrec.xi2 = dscomResult.xi2;
                satrec.xi3 = dscomResult.xi3;
                satrec.xl2 = dscomResult.xl2;
                satrec.xl3 = dscomResult.xl3;
                satrec.xl4 = dscomResult.xl4;
                satrec.zmol = dscomResult.zmol;
                satrec.zmos = dscomResult.zmos;
                const { nm, z1, z3, z11, z13, z21, z23, z31, z33 } = dscomResult;
                const dpperOptions = {
                    inclo: inclm,
                    init: satrec.init,
                    ep: satrec.ecco,
                    inclp: satrec.inclo,
                    nodep: satrec.nodeo,
                    argpp: satrec.argpo,
                    mp: satrec.mo,
                    opsmode: satrec.operationmode,
                    satrec,
                };
                const dpperResult = Sgp4.dpper_(dpperOptions);
                satrec.ecco = dpperResult.ep;
                satrec.inclo = dpperResult.inclp;
                satrec.nodeo = dpperResult.nodep;
                satrec.argpo = dpperResult.argpp;
                satrec.mo = dpperResult.mp;
                const argpm = 0.0;
                const nodem = 0.0;
                const mm = 0.0;
                const dsinitOptions = {
                    xke,
                    cosim,
                    emsq,
                    argpo: satrec.argpo,
                    s1,
                    s2,
                    s3,
                    s4,
                    s5,
                    sinim,
                    ss1,
                    ss2,
                    ss3,
                    ss4,
                    ss5,
                    sz1,
                    sz3,
                    sz11,
                    sz13,
                    sz21,
                    sz23,
                    sz31,
                    sz33,
                    t: satrec.t,
                    tc,
                    gsto: satrec.gsto,
                    mo: satrec.mo,
                    mdot: satrec.mdot,
                    no: satrec.no,
                    nodeo: satrec.nodeo,
                    nodedot: satrec.nodedot,
                    xPIdot,
                    z1,
                    z3,
                    z11,
                    z13,
                    z21,
                    z23,
                    z31,
                    z33,
                    ecco: satrec.ecco,
                    eccsq,
                    em,
                    argpm,
                    inclm,
                    mm,
                    nm,
                    nodem,
                    irez: satrec.irez,
                    atime: satrec.atime,
                    d2201: satrec.d2201,
                    d2211: satrec.d2211,
                    d3210: satrec.d3210,
                    d3222: satrec.d3222,
                    d4410: satrec.d4410,
                    d4422: satrec.d4422,
                    d5220: satrec.d5220,
                    d5232: satrec.d5232,
                    d5421: satrec.d5421,
                    d5433: satrec.d5433,
                    dedt: satrec.dedt,
                    didt: satrec.didt,
                    dmdt: satrec.dmdt,
                    dnodt: satrec.dnodt,
                    domdt: satrec.domdt,
                    del1: satrec.del1,
                    del2: satrec.del2,
                    del3: satrec.del3,
                    xfact: satrec.xfact,
                    xlamo: satrec.xlamo,
                    xli: satrec.xli,
                    xni: satrec.xni,
                };
                const dsinitResult = Sgp4.dsinit_(dsinitOptions);
                satrec.irez = dsinitResult.irez;
                satrec.atime = dsinitResult.atime;
                satrec.d2201 = dsinitResult.d2201;
                satrec.d2211 = dsinitResult.d2211;
                satrec.d3210 = dsinitResult.d3210;
                satrec.d3222 = dsinitResult.d3222;
                satrec.d4410 = dsinitResult.d4410;
                satrec.d4422 = dsinitResult.d4422;
                satrec.d5220 = dsinitResult.d5220;
                satrec.d5232 = dsinitResult.d5232;
                satrec.d5421 = dsinitResult.d5421;
                satrec.d5433 = dsinitResult.d5433;
                satrec.dedt = dsinitResult.dedt;
                satrec.didt = dsinitResult.didt;
                satrec.dmdt = dsinitResult.dmdt;
                satrec.dnodt = dsinitResult.dnodt;
                satrec.domdt = dsinitResult.domdt;
                satrec.del1 = dsinitResult.del1;
                satrec.del2 = dsinitResult.del2;
                satrec.del3 = dsinitResult.del3;
                satrec.xfact = dsinitResult.xfact;
                satrec.xlamo = dsinitResult.xlamo;
                satrec.xli = dsinitResult.xli;
                satrec.xni = dsinitResult.xni;
            }
            // ----------- set variables if not deep space -----------
            if (!satrec.isimp) {
                const cc1sq = satrec.cc1 * satrec.cc1;
                satrec.d2 = 4.0 * ao * tsi * cc1sq;
                const temp = (satrec.d2 * tsi * satrec.cc1) / 3.0;
                satrec.d3 = (17.0 * ao + sfour) * temp;
                satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * satrec.cc1;
                satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
                satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 * (12.0 * satrec.d2 + 10.0 * cc1sq));
                satrec.t5cof =
                    0.2 *
                        (3.0 * satrec.d4 +
                            12.0 * satrec.cc1 * satrec.d3 +
                            6.0 * satrec.d2 * satrec.d2 +
                            15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
            } // If omeosq = 0 ...
            /* Finally propogate to zero epoch to initialize all others. */
            /*
             * Sgp4fix take out check to let satellites process until they are actually below earth surface
             * if(satrec.error == 0)
             */
        }
        Sgp4.propagate(satrec, 0);
        satrec.init = false;
        /*
         * Sgp4fix return boolean. satrec.error contains any error codes
         *  return satrec; -- no reason to return anything in JS
         */
    }
}


}),
"./src/engine/ootk/src/time/Epoch.ts": 
/*!*******************************************!*\
  !*** ./src/engine/ootk/src/time/Epoch.ts ***!
  \*******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Epoch: () => (Epoch)
});
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

// / Base class for [Epoch] data.
class Epoch {
    posix;
    /*
     * Create a new [Epoch] object given the number of seconds elapsed since the
     * [posix] epoch _(`1970-01-01T00:00:00.000`)_ in the [Epoch] time scale.
     */
    constructor(posix) {
        this.posix = posix;
        if (posix < 0) {
            throw new Error('Epoch cannot be negative');
        }
    }
    toString() {
        return this.toDateTime().toISOString();
    }
    // / Convert this to an Excel spreadsheet string.
    toExcelString() {
        return this.toString().substring(0, 19);
    }
    // / Return the difference _(s)_ between this and another [epoch]/
    difference(epoch) {
        return this.posix - epoch.posix;
    }
    // / Check if this has the same timestamp as the provided [epoch].
    equals(epoch) {
        return this.posix === epoch.posix;
    }
    // / Convert to a [DateTime] object.
    toDateTime() {
        return new Date(this.posix * 1000);
    }
    toEpochYearAndDay() {
        const currentDateObj = this.toDateTime();
        const epochYear = currentDateObj.getUTCFullYear().toString().slice(2, 4);
        const epochDay = this.getDayOfYear_(currentDateObj);
        const timeOfDay = (currentDateObj.getUTCHours() * 60 + currentDateObj.getUTCMinutes()) / 1440;
        const epochDayStr = (epochDay + timeOfDay).toFixed(8).padStart(12, '0');
        return {
            epochYr: epochYear,
            epochDay: epochDayStr,
        };
    }
    getDayOfYear_(date) {
        const dayCount = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
        const mn = date.getUTCMonth();
        const dn = date.getUTCDate();
        let dayOfYear = dayCount[mn] + dn;
        if (mn > 1 && this.isLeapYear_(date)) {
            dayOfYear++;
        }
        return dayOfYear;
    }
    isLeapYear_(dateIn) {
        const year = dateIn.getUTCFullYear();
        if ((year & 3) !== 0) {
            return false;
        }
        return year % 100 !== 0 || year % 400 === 0;
    }
    // / Convert to Julian date.
    toJulianDate() {
        return this.posix / _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerDay + 2440587.5;
    }
    // / Convert to Julian centuries.
    toJulianCenturies() {
        return (this.toJulianDate() - 2451545) / 36525;
    }
    // / Check if this is later than the [other] epoch.
    operatorGreaterThan(other) {
        return this.posix > other.posix;
    }
    // / Check if this is later or the same as the [other] epoch.
    operatorGreaterThanOrEqual(other) {
        return this.posix >= other.posix;
    }
    // / Check if this is earlier than the [other] epoch.
    operatorLessThan(other) {
        return this.posix < other.posix;
    }
    // / Check if this is earlier or the same as the [other] epoch.
    operatorLessThanOrEqual(other) {
        return this.posix <= other.posix;
    }
}


}),
"./src/engine/ootk/src/time/EpochGPS.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/time/EpochGPS.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochGPS: () => (EpochGPS)
});
/* ESM import */var _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../data/DataHandler.js */ "./src/engine/ootk/src/data/DataHandler.ts");
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */


// / Global Positioning System _(GPS)_ formatted epoch.
class EpochGPS {
    week;
    seconds;
    /**
     * Create a new GPS epoch given the [week] since reference epoch, and number
     * of [seconds] into the [week].
     * @param week Number of weeks since the GPS reference epoch.
     * @param seconds Number of seconds into the week.
     * @param reference Reference should always be EpochUTC.fromDateTimeString('1980-01-06T00:00:00.000Z').
     */
    constructor(week, seconds, reference) {
        this.week = week;
        this.seconds = seconds;
        if (week < 0) {
            throw new Error('GPS week must be non-negative.');
        }
        if (seconds < 0 || seconds >= _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.secondsPerWeek) {
            throw new Error('GPS seconds must be within a week.');
        }
        // TODO: Set EpochGPS.reference statically without circular dependency.
        EpochGPS.reference = reference;
    }
    // / Number of weeks since the GPS reference epoch.
    static reference;
    // / GPS leap second difference from TAI/UTC offsets.
    static offset = 19;
    // / Get GPS week accounting for 10-bit rollover.
    get week10Bit() {
        return this.week % 2 ** 10;
    }
    // / Get GPS week accounting for 13-bit rollover.
    get week13Bit() {
        return this.week % 2 ** 13;
    }
    toString() {
        return `${this.week}:${this.seconds.toFixed(3)}`;
    }
    // / Convert this to a UTC epoch.
    toUTC() {
        const init = EpochGPS.reference.roll((this.week * _utils_constants_js__WEBPACK_IMPORTED_MODULE_1__.secondsPerWeek + this.seconds));
        const ls = _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_0__.DataHandler.getInstance().getLeapSeconds(init.toJulianDate());
        return init.roll(-(ls - EpochGPS.offset));
    }
}


}),
"./src/engine/ootk/src/time/EpochTAI.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/time/EpochTAI.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochTAI: () => (EpochTAI)
});
/* ESM import */var _Epoch_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Epoch.js */ "./src/engine/ootk/src/time/Epoch.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/** Represents an Epoch in International Atomic Time (TAI). */
class EpochTAI extends _Epoch_js__WEBPACK_IMPORTED_MODULE_0__.Epoch {
}


}),
"./src/engine/ootk/src/time/EpochTDB.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/time/EpochTDB.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochTDB: () => (EpochTDB)
});
/* ESM import */var _Epoch_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Epoch.js */ "./src/engine/ootk/src/time/Epoch.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/** Represents an Epoch in TDB (Barycentric Dynamical Time). */
class EpochTDB extends _Epoch_js__WEBPACK_IMPORTED_MODULE_0__.Epoch {
}


}),
"./src/engine/ootk/src/time/EpochTT.ts": 
/*!*********************************************!*\
  !*** ./src/engine/ootk/src/time/EpochTT.ts ***!
  \*********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochTT: () => (EpochTT)
});
/* ESM import */var _Epoch_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Epoch.js */ "./src/engine/ootk/src/time/Epoch.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/** Represents a Terrestrial Time (TT) epoch. */
class EpochTT extends _Epoch_js__WEBPACK_IMPORTED_MODULE_0__.Epoch {
}


}),
"./src/engine/ootk/src/time/EpochUTC.ts": 
/*!**********************************************!*\
  !*** ./src/engine/ootk/src/time/EpochUTC.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochUTC: () => (EpochUTC)
});
/* ESM import */var _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../utils/constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _utils_functions_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../utils/functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./../data/DataHandler.js */ "./src/engine/ootk/src/data/DataHandler.ts");
/* ESM import */var _Epoch_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./Epoch.js */ "./src/engine/ootk/src/time/Epoch.ts");
/* ESM import */var _EpochGPS_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./EpochGPS.js */ "./src/engine/ootk/src/time/EpochGPS.ts");
/* ESM import */var _EpochTAI_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./EpochTAI.js */ "./src/engine/ootk/src/time/EpochTAI.ts");
/* ESM import */var _EpochTDB_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./EpochTDB.js */ "./src/engine/ootk/src/time/EpochTDB.ts");
/* ESM import */var _EpochTT_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./EpochTT.js */ "./src/engine/ootk/src/time/EpochTT.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */








class EpochUTC extends _Epoch_js__WEBPACK_IMPORTED_MODULE_3__.Epoch {
    static now() {
        return new EpochUTC(new Date().getTime() / 1000);
    }
    static fromDate({ year, month, day, hour = 0, minute = 0, second = 0 }) {
        return new EpochUTC(EpochUTC.dateToPosix_({ year, month, day, hour, minute, second }));
    }
    static fromDateTime(dt) {
        return new EpochUTC(dt.getTime() / 1000);
    }
    static fromDateTimeString(dateTimeString) {
        const dts = dateTimeString.trim().toUpperCase().endsWith('Z') ? dateTimeString : `${dateTimeString}Z`;
        return new EpochUTC(new Date(dts).getTime() / 1000);
    }
    static fromJ2000TTSeconds(seconds) {
        const tInit = new EpochUTC(seconds + 946728000);
        const ls = _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_2__.DataHandler.getInstance().getLeapSeconds(tInit.toJulianDate());
        return tInit.roll(-32.184 - ls);
    }
    static fromDefinitiveString(definitiveString) {
        const fields = definitiveString.trim().split(' ');
        const dateFields = fields[0].split('/');
        const day = parseInt(dateFields[0]);
        const year = parseInt(dateFields[1]);
        // eslint-disable-next-line prefer-destructuring
        const timeField = fields[1];
        // Add day - 1 days in milliseconds to the epoch.
        const dts = new Date(`${year}-01-01T${timeField}Z`).getTime() + (day - 1) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.MS_PER_DAY;
        return new EpochUTC(dts / 1000);
    }
    roll(seconds) {
        return new EpochUTC(this.posix + seconds);
    }
    toMjd() {
        return this.toJulianDate() - 2400000.5;
    }
    toMjdGsfc() {
        return this.toMjd() - 29999.5;
    }
    toTAI() {
        const ls = _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_2__.DataHandler.getInstance().getLeapSeconds(this.toJulianDate());
        return new _EpochTAI_js__WEBPACK_IMPORTED_MODULE_5__.EpochTAI(this.posix + ls);
    }
    toTT() {
        return new _EpochTT_js__WEBPACK_IMPORTED_MODULE_7__.EpochTT(this.toTAI().posix + 32.184);
    }
    toTDB() {
        const tt = this.toTT();
        const tTT = tt.toJulianCenturies();
        const mEarth = (357.5277233 + 35999.05034 * tTT) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
        const seconds = 0.001658 * Math.sin(mEarth) + 0.00001385 * Math.sin(2 * mEarth);
        return new _EpochTDB_js__WEBPACK_IMPORTED_MODULE_6__.EpochTDB(tt.posix + seconds);
    }
    toGPS() {
        const referenceTime = EpochUTC.fromDateTimeString('1980-01-06T00:00:00.000Z');
        const ls = _data_DataHandler_js__WEBPACK_IMPORTED_MODULE_2__.DataHandler.getInstance().getLeapSeconds(this.toJulianDate());
        const delta = this.roll(ls - _EpochGPS_js__WEBPACK_IMPORTED_MODULE_4__.EpochGPS.offset).difference(referenceTime);
        const week = delta / _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerWeek;
        const weekFloor = Math.floor(week);
        const seconds = (week - weekFloor) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.secondsPerWeek;
        return new _EpochGPS_js__WEBPACK_IMPORTED_MODULE_4__.EpochGPS(weekFloor, seconds, referenceTime);
    }
    gmstAngle() {
        const t = this.toJulianCenturies();
        const seconds = (0,_utils_functions_js__WEBPACK_IMPORTED_MODULE_1__.evalPoly)(t, EpochUTC.gmstPoly_);
        let result = ((seconds / 240) * _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD) % _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
        if (result < 0) {
            result += _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
        }
        return result;
    }
    gmstAngleDegrees() {
        return this.gmstAngle() * _utils_constants_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG;
    }
    static gmstPoly_ = new Float64Array([
        -6.2e-6,
        0.093104,
        876600 * 3600 + 8640184.812866,
        67310.54841,
    ]);
    static dayOfYearLookup_ = [
        [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
        [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335],
    ];
    static isLeapYear_(year) {
        return (year % 4 === 0 && year % 100 !== 0) || year % 400 === 0;
    }
    static dayOfYear_(year, month, day) {
        const dex = EpochUTC.isLeapYear_(year) ? 1 : 0;
        const dayOfYearArray = EpochUTC.dayOfYearLookup_[dex];
        const daysBeforeCurrentMonth = dayOfYearArray[month - 1];
        return daysBeforeCurrentMonth + day - 1;
    }
    static dateToPosix_({ year, month, day, hour, minute, second }) {
        const days = EpochUTC.dayOfYear_(year, month, day);
        const yearMod = year - 1900;
        return (minute * 60 +
            hour * 3600 +
            days * 86400 +
            (yearMod - 70) * 31536000 +
            Math.floor((yearMod - 69) / 4) * 86400 -
            Math.floor((yearMod - 1) / 100) * 86400 +
            Math.floor((yearMod + 299) / 400) * 86400 +
            second);
    }
}


}),
"./src/engine/ootk/src/time/EpochWindow.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/time/EpochWindow.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  EpochWindow: () => (EpochWindow)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
class EpochWindow {
    start;
    end;
    constructor(start, end) {
        this.start = start;
        this.end = end;
    }
}


}),
"./src/engine/ootk/src/time/TimeStamped.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/time/TimeStamped.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  TimeStamped: () => (TimeStamped)
});
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
class TimeStamped {
    /**
     * Create a new time stamped value container at the provided epoch.
     * @param epoch The timestamp epoch.
     * @param value The timestamped value.
     */
    constructor(epoch, value) {
        this.epoch = epoch;
        this.value = value;
    }
    /**
     * Timestamp epoch.
     */
    epoch;
    /**
     * Timestamped value.
     */
    value;
}


}),
"./src/engine/ootk/src/time/index.ts": 
/*!*******************************************!*\
  !*** ./src/engine/ootk/src/time/index.ts ***!
  \*******************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  Epoch: () => (/* reexport safe */ _Epoch_js__WEBPACK_IMPORTED_MODULE_0__.Epoch),
  EpochGPS: () => (/* reexport safe */ _EpochGPS_js__WEBPACK_IMPORTED_MODULE_1__.EpochGPS),
  EpochTAI: () => (/* reexport safe */ _EpochTAI_js__WEBPACK_IMPORTED_MODULE_2__.EpochTAI),
  EpochTDB: () => (/* reexport safe */ _EpochTDB_js__WEBPACK_IMPORTED_MODULE_3__.EpochTDB),
  EpochTT: () => (/* reexport safe */ _EpochTT_js__WEBPACK_IMPORTED_MODULE_4__.EpochTT),
  EpochUTC: () => (/* reexport safe */ _EpochUTC_js__WEBPACK_IMPORTED_MODULE_5__.EpochUTC),
  EpochWindow: () => (/* reexport safe */ _EpochWindow_js__WEBPACK_IMPORTED_MODULE_6__.EpochWindow),
  TimeStamped: () => (/* reexport safe */ _TimeStamped_js__WEBPACK_IMPORTED_MODULE_7__.TimeStamped)
});
/* ESM import */var _Epoch_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./Epoch.js */ "./src/engine/ootk/src/time/Epoch.ts");
/* ESM import */var _EpochGPS_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./EpochGPS.js */ "./src/engine/ootk/src/time/EpochGPS.ts");
/* ESM import */var _EpochTAI_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./EpochTAI.js */ "./src/engine/ootk/src/time/EpochTAI.ts");
/* ESM import */var _EpochTDB_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./EpochTDB.js */ "./src/engine/ootk/src/time/EpochTDB.ts");
/* ESM import */var _EpochTT_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./EpochTT.js */ "./src/engine/ootk/src/time/EpochTT.ts");
/* ESM import */var _EpochUTC_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(/*! ./EpochUTC.js */ "./src/engine/ootk/src/time/EpochUTC.ts");
/* ESM import */var _EpochWindow_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(/*! ./EpochWindow.js */ "./src/engine/ootk/src/time/EpochWindow.ts");
/* ESM import */var _TimeStamped_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(/*! ./TimeStamped.js */ "./src/engine/ootk/src/time/TimeStamped.ts");










}),
"./src/engine/ootk/src/transforms/conversions.ts": 
/*!*******************************************************!*\
  !*** ./src/engine/ootk/src/transforms/conversions.ts ***!
  \*******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  deg2rad: () => (deg2rad),
  getDegLat: () => (getDegLat),
  getDegLon: () => (getDegLon),
  getRadLat: () => (getRadLat),
  getRadLon: () => (getRadLon),
  rad2deg: () => (rad2deg)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");
/**
 * @author Theodore Kruczek
 * @description Orbital Object ToolKit (ootk) is a collection of tools for working
 * with satellites and other orbital objects.
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Many of the classes are based off of the work of @david-rc-dayton and his
 * Pious Squid library (https://github.com/david-rc-dayton/pious_squid) which
 * is licensed under the MIT license.
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Converts radians to degrees.
 * @param radians The value in radians to be converted.
 * @returns The value in degrees.
 */
function rad2deg(radians) {
    return ((radians * 180) / _main_js__WEBPACK_IMPORTED_MODULE_0__.PI);
}
/**
 * Converts degrees to radians.
 * @param degrees The value in degrees to be converted.
 * @returns The value in radians.
 */
function deg2rad(degrees) {
    return ((degrees * _main_js__WEBPACK_IMPORTED_MODULE_0__.PI) / 180.0);
}
/**
 * Converts radians to degrees latitude.
 * @param radians The radians value to convert.
 * @returns The corresponding degrees latitude.
 * @throws RangeError if the radians value is outside the range [-PI/2; PI/2].
 */
function getDegLat(radians) {
    if (radians < -_main_js__WEBPACK_IMPORTED_MODULE_0__.PI / 2 || radians > _main_js__WEBPACK_IMPORTED_MODULE_0__.PI / 2) {
        throw new RangeError('Latitude radians must be in range [-PI/2; PI/2].');
    }
    return (radians * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
}
/**
 * Converts radians to degrees for longitude.
 * @param radians The value in radians to be converted.
 * @returns The converted value in degrees.
 * @throws {RangeError} If the input radians is not within the range [-PI; PI].
 */
function getDegLon(radians) {
    if (radians < -_main_js__WEBPACK_IMPORTED_MODULE_0__.PI || radians > _main_js__WEBPACK_IMPORTED_MODULE_0__.PI) {
        throw new RangeError('Longitude radians must be in range [-PI; PI].');
    }
    return (radians * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
}
/**
 * Converts degrees to radians for latitude.
 * @param degrees The degrees value to convert.
 * @returns The equivalent radians value.
 * @throws {RangeError} If the degrees value is not within the range [-90, 90].
 */
function getRadLat(degrees) {
    if (degrees < -90 || degrees > 90) {
        throw new RangeError('Latitude degrees must be in range [-90; 90].');
    }
    return (degrees * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
}
/**
 * Converts degrees to radians.
 * @param degrees The value in degrees to be converted.
 * @returns The value in radians.
 * @throws {RangeError} If the input degrees are not within the range [-180; 180].
 */
function getRadLon(degrees) {
    if (degrees < -180 || degrees > 180) {
        throw new RangeError('Longitude degrees must be in range [-180; 180].');
    }
    return (degrees * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
}


}),
"./src/engine/ootk/src/transforms/index.ts": 
/*!*************************************************!*\
  !*** ./src/engine/ootk/src/transforms/index.ts ***!
  \*************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  azel2uv: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.azel2uv),
  calcGmst: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.calcGmst),
  calcIncFromAz: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.calcIncFromAz),
  calcInertAz: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.calcInertAz),
  deg2rad: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.deg2rad),
  ecf2eci: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.ecf2eci),
  ecf2enu: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.ecf2enu),
  ecf2rae: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.ecf2rae),
  ecfRad2rae: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.ecfRad2rae),
  eci2ecf: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.eci2ecf),
  eci2lla: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.eci2lla),
  eci2rae: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.eci2rae),
  enu2rf: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.enu2rf),
  getDegLat: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.getDegLat),
  getDegLon: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.getDegLon),
  getRadLat: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.getRadLat),
  getRadLon: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.getRadLon),
  jday: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.jday),
  lla2ecef: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.lla2ecef),
  lla2ecf: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.lla2ecf),
  lla2eci: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.lla2eci),
  lla2sez: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.lla2sez),
  llaRad2ecf: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.llaRad2ecf),
  rad2deg: () => (/* reexport safe */ _conversions_js__WEBPACK_IMPORTED_MODULE_0__.rad2deg),
  rae2ecf: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2ecf),
  rae2eci: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2eci),
  rae2enu: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2enu),
  rae2raeOffBoresight: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2raeOffBoresight),
  rae2ruv: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2ruv),
  rae2sez: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.rae2sez),
  sez2rae: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.sez2rae),
  uv2azel: () => (/* reexport safe */ _transforms_js__WEBPACK_IMPORTED_MODULE_1__.uv2azel)
});
/* ESM import */var _conversions_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./conversions.js */ "./src/engine/ootk/src/transforms/conversions.ts");
/* ESM import */var _transforms_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./transforms.js */ "./src/engine/ootk/src/transforms/transforms.ts");




}),
"./src/engine/ootk/src/transforms/transforms.ts": 
/*!******************************************************!*\
  !*** ./src/engine/ootk/src/transforms/transforms.ts ***!
  \******************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  azel2uv: () => (azel2uv),
  calcGmst: () => (calcGmst),
  calcIncFromAz: () => (calcIncFromAz),
  calcInertAz: () => (calcInertAz),
  ecf2eci: () => (ecf2eci),
  ecf2enu: () => (ecf2enu),
  ecf2rae: () => (ecf2rae),
  ecfRad2rae: () => (ecfRad2rae),
  eci2ecf: () => (eci2ecf),
  eci2lla: () => (eci2lla),
  eci2rae: () => (eci2rae),
  enu2rf: () => (enu2rf),
  jday: () => (jday),
  lla2ecef: () => (lla2ecef),
  lla2ecf: () => (lla2ecf),
  lla2eci: () => (lla2eci),
  lla2sez: () => (lla2sez),
  llaRad2ecf: () => (llaRad2ecf),
  rae2ecf: () => (rae2ecf),
  rae2eci: () => (rae2eci),
  rae2enu: () => (rae2enu),
  rae2raeOffBoresight: () => (rae2raeOffBoresight),
  rae2ruv: () => (rae2ruv),
  rae2sez: () => (rae2sez),
  sez2rae: () => (sez2rae),
  uv2azel: () => (uv2azel)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");

/**
 * Converts ECF to ECI coordinates.
 *
 * [X]     [C -S  0][X]
 * [Y]  =  [S  C  0][Y]
 * [Z]eci  [0  0  1][Z]ecf
 * @param ecf takes xyz coordinates
 * @param gmst takes a number in gmst time
 * @returns array containing eci coordinates
 */
function ecf2eci(ecf, gmst) {
    const X = (ecf.x * Math.cos(gmst) - ecf.y * Math.sin(gmst));
    const Y = (ecf.x * Math.sin(gmst) + ecf.y * Math.cos(gmst));
    const Z = ecf.z;
    return { x: X, y: Y, z: Z };
}
/**
 * Converts ECEF coordinates to ENU coordinates.
 * @param ecf - The ECEF coordinates.
 * @param lla - The LLA coordinates.
 * @returns The ENU coordinates.
 */
function ecf2enu(ecf, lla) {
    const { lat, lon } = lla;
    const { x, y, z } = ecf;
    const e = (-Math.sin(lon) * x + Math.cos(lon) * y);
    const n = (-Math.sin(lat) * Math.cos(lon) * x - Math.sin(lat) * Math.sin(lon) * y + Math.cos(lat) * z);
    const u = (Math.cos(lat) * Math.cos(lon) * x + Math.cos(lat) * Math.sin(lon) * y + Math.sin(lat) * z);
    return { x: e, y: n, z: u };
}
/**
 * Converts ECI to ECF coordinates.
 *
 * [X]     [C -S  0][X]
 * [Y]  =  [S  C  0][Y]
 * [Z]eci  [0  0  1][Z]ecf
 *
 * Inverse:
 * [X]     [C  S  0][X]
 * [Y]  =  [-S C  0][Y]
 * [Z]ecf  [0  0  1][Z]eci
 * @param eci takes xyz coordinates
 * @param gmst takes a number in gmst time
 * @returns array containing ecf coordinates
 */
function eci2ecf(eci, gmst) {
    const x = (eci.x * Math.cos(gmst) + eci.y * Math.sin(gmst));
    const y = (eci.x * -Math.sin(gmst) + eci.y * Math.cos(gmst));
    const z = eci.z;
    return {
        x,
        y,
        z,
    };
}
/**
 * EciToGeodetic converts eci coordinates to lla coordinates
 * @variation cached - results are cached
 * @param eci takes xyz coordinates
 * @param gmst takes a number in gmst time
 * @returns array containing lla coordinates
 */
function eci2lla(eci, gmst) {
    // http://www.celestrak.com/columns/v02n03/
    const a = 6378.137;
    const b = 6356.7523142;
    const R = Math.sqrt(eci.x * eci.x + eci.y * eci.y);
    const f = (a - b) / a;
    const e2 = 2 * f - f * f;
    let lon = Math.atan2(eci.y, eci.x) - gmst;
    while (lon < -_main_js__WEBPACK_IMPORTED_MODULE_0__.PI) {
        lon += _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
    }
    while (lon > _main_js__WEBPACK_IMPORTED_MODULE_0__.PI) {
        lon -= _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU;
    }
    const kmax = 20;
    let k = 0;
    let lat = Math.atan2(eci.z, Math.sqrt(eci.x * eci.x + eci.y * eci.y));
    let C = 0;
    while (k < kmax) {
        C = 1 / Math.sqrt(1 - e2 * (Math.sin(lat) * Math.sin(lat)));
        lat = Math.atan2(eci.z + a * C * e2 * Math.sin(lat), R);
        k += 1;
    }
    const alt = R / Math.cos(lat) - a * C;
    lon = (lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
    lat = (lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
    return { lon: lon, lat: lat, alt: alt };
}
/**
 * Converts geodetic coordinates (longitude, latitude, altitude) to Earth-Centered Earth-Fixed (ECF) coordinates.
 * @param lla The geodetic coordinates in radians and meters.
 * @returns The ECF coordinates in meters.
 */
function llaRad2ecf(lla) {
    const { lon, lat, alt } = lla;
    const a = 6378.137;
    const b = 6356.7523142;
    const f = (a - b) / a;
    const e2 = 2 * f - f * f;
    const normal = a / Math.sqrt(1 - e2 * Math.sin(lat) ** 2);
    const x = (normal + alt) * Math.cos(lat) * Math.cos(lon);
    const y = (normal + alt) * Math.cos(lat) * Math.sin(lon);
    const z = (normal * (1 - e2) + alt) * Math.sin(lat);
    return {
        x: x,
        y: y,
        z: z,
    };
}
/**
 * Converts geodetic coordinates (longitude, latitude, altitude) to Earth-Centered Earth-Fixed (ECF) coordinates.
 * @param lla The geodetic coordinates in degrees and meters.
 * @returns The ECF coordinates in meters.
 */
function lla2ecf(lla) {
    const { lon, lat, alt } = lla;
    const lonRad = lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    const latRad = lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    return llaRad2ecf({
        lon: lonRad,
        lat: latRad,
        alt,
    });
}
/**
 * Converts geodetic coordinates (latitude, longitude, altitude) to Earth-centered inertial (ECI) coordinates.
 * @variation cached - results are cached
 * @param lla The geodetic coordinates in radians and meters.
 * @param gmst The Greenwich Mean Sidereal Time in seconds.
 * @returns The ECI coordinates in meters.
 */
function lla2eci(lla, gmst) {
    const { lat, lon, alt } = lla;
    const cosLat = Math.cos(lat);
    const sinLat = Math.sin(lat);
    const cosLon = Math.cos(lon + gmst);
    const sinLon = Math.sin(lon + gmst);
    const x = (_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean + alt) * cosLat * cosLon;
    const y = (_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean + alt) * cosLat * sinLon;
    const z = (_main_js__WEBPACK_IMPORTED_MODULE_0__.Earth.radiusMean + alt) * sinLat;
    return { x, y, z };
}
/**
 * Calculates Geodetic Lat Lon Alt to ECEF coordinates.
 * @deprecated This needs to be validated.
 * @param lla The geodetic coordinates in degrees and meters.
 * @returns The ECEF coordinates in meters.
 */
function lla2ecef(lla) {
    const { lat, lon, alt } = lla;
    const a = 6378.137; // semi-major axis length in meters according to the WGS84
    const b = 6356.752314245; // semi-minor axis length in meters according to the WGS84
    const e = Math.sqrt(1 - b ** 2 / a ** 2); // eccentricity
    const N = a / Math.sqrt(1 - e ** 2 * Math.sin(lat) ** 2); // radius of curvature in the prime vertical
    const x = ((N + alt) * Math.cos(lat) * Math.cos(lon));
    const y = ((N + alt) * Math.cos(lat) * Math.sin(lon));
    const z = ((N * (1 - e ** 2) + alt) * Math.sin(lat));
    return { x, y, z };
}
/**
 * Converts LLA to SEZ coordinates.
 * @see http://www.celestrak.com/columns/v02n02/
 * @param lla The LLA coordinates.
 * @param ecf The ECF coordinates.
 * @returns The SEZ coordinates.
 */
function lla2sez(lla, ecf) {
    const lon = lla.lon;
    const lat = lla.lat;
    const observerEcf = llaRad2ecf({
        lat,
        lon,
        alt: 0,
    });
    const rx = ecf.x - observerEcf.x;
    const ry = ecf.y - observerEcf.y;
    const rz = ecf.z - observerEcf.z;
    // Top is short for topocentric
    const south = Math.sin(lat) * Math.cos(lon) * rx + Math.sin(lat) * Math.sin(lon) * ry - Math.cos(lat) * rz;
    const east = -Math.sin(lon) * rx + Math.cos(lon) * ry;
    const zenith = Math.cos(lat) * Math.cos(lon) * rx + Math.cos(lat) * Math.sin(lon) * ry + Math.sin(lat) * rz;
    return { s: south, e: east, z: zenith };
}
/**
 * Converts a vector in Right Ascension, Elevation, and Range (RAE) coordinate system
 * to a vector in South, East, and Zenith (SEZ) coordinate system.
 * @param rae The vector in RAE coordinate system.
 * @returns The vector in SEZ coordinate system.
 */
function rae2sez(rae) {
    const south = -rae.rng * Math.cos(rae.el) * Math.cos(rae.az);
    const east = rae.rng * Math.cos(rae.el) * Math.sin(rae.az);
    const zenith = rae.rng * Math.sin(rae.el);
    return {
        s: south,
        e: east,
        z: zenith,
    };
}
/**
 * Converts a vector in Right Ascension, Elevation, and Range (RAE) coordinate system
 * to Earth-Centered Fixed (ECF) coordinate system.
 * @template D - The dimension of the RAE vector.
 * @template A - The dimension of the LLA vector.
 * @param rae - The vector in RAE coordinate system.
 * @param lla - The vector in LLA coordinate system.
 * @returns The vector in ECF coordinate system.
 */
function rae2ecf(rae, lla) {
    const llaRad = {
        lat: (lla.lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        lon: (lla.lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        alt: lla.alt,
    };
    const raeRad = {
        az: (rae.az * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        el: (rae.el * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        rng: rae.rng,
    };
    const obsEcf = llaRad2ecf(llaRad);
    const sez = rae2sez(raeRad);
    // Some needed calculations
    const slat = Math.sin(llaRad.lat);
    const slon = Math.sin(llaRad.lon);
    const clat = Math.cos(llaRad.lat);
    const clon = Math.cos(llaRad.lon);
    const x = slat * clon * sez.s + -slon * sez.e + clat * clon * sez.z + obsEcf.x;
    const y = slat * slon * sez.s + clon * sez.e + clat * slon * sez.z + obsEcf.y;
    const z = -clat * sez.s + slat * sez.z + obsEcf.z;
    return { x, y, z };
}
/**
 * Converts a vector from RAE (Range, Azimuth, Elevation) coordinates to ECI (Earth-Centered Inertial) coordinates.
 * @variation cached - results are cached
 * @param rae The vector in RAE coordinates.
 * @param lla The vector in LLA (Latitude, Longitude, Altitude) coordinates.
 * @param gmst The Greenwich Mean Sidereal Time.
 * @returns The vector in ECI coordinates.
 */
function rae2eci(rae, lla, gmst) {
    const ecf = rae2ecf(rae, lla);
    const eci = ecf2eci(ecf, gmst);
    return eci;
}
/**
 * Converts a vector in RAE (Range, Azimuth, Elevation) coordinates to ENU (East, North, Up) coordinates.
 * @param rae - The vector in RAE coordinates.
 * @returns The vector in ENU coordinates.
 */
function rae2enu(rae) {
    const e = (rae.rng * Math.cos(rae.el) * Math.sin(rae.az));
    const n = (rae.rng * Math.cos(rae.el) * Math.cos(rae.az));
    const u = (rae.rng * Math.sin(rae.el));
    return { x: e, y: n, z: u };
}
/**
 * Converts South, East, and Zenith (SEZ) coordinates to Right Ascension, Elevation, and Range (RAE) coordinates.
 * @param sez The SEZ coordinates.
 * @returns Rng, Az, El array
 */
function sez2rae(sez) {
    const rng = Math.sqrt(sez.s * sez.s + sez.e * sez.e + sez.z * sez.z);
    const el = Math.asin(sez.z / rng);
    const az = (Math.atan2(-sez.e, sez.s) + _main_js__WEBPACK_IMPORTED_MODULE_0__.PI);
    return { rng, az, el };
}
/**
 * Converts Earth-Centered Fixed (ECF) coordinates to Right Ascension (RA),
 * Elevation (E), and Azimuth (A) coordinates.
 * @param lla The Latitude, Longitude, and Altitude (LLA) coordinates.
 * @param ecf The Earth-Centered Fixed (ECF) coordinates.
 * @returns The Right Ascension (RA), Elevation (E), and Azimuth (A) coordinates.
 */
function ecfRad2rae(lla, ecf) {
    const sezCoords = lla2sez(lla, ecf);
    const rae = sez2rae(sezCoords);
    return { rng: rae.rng, az: (rae.az * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG), el: (rae.el * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG) };
}
/**
 * Converts Earth-Centered Fixed (ECF) coordinates to Right Ascension (RA),
 * Elevation (E), and Azimuth (A) coordinates.
 * @variation cached - results are cached
 * @param lla The Latitude, Longitude, and Altitude (LLA) coordinates.
 * @param ecf The Earth-Centered Fixed (ECF) coordinates.
 * @returns The Right Ascension (RA), Elevation (E), and Azimuth (A) coordinates.
 */
function ecf2rae(lla, ecf) {
    const { lat, lon } = lla;
    const latRad = (lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    const lonRad = (lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    const rae = ecfRad2rae({ lat: latRad, lon: lonRad, alt: lla.alt }, ecf);
    return rae;
}
const jday = (year, mon, day, hr, minute, sec) => {
    if (typeof year === 'undefined') {
        const now = new Date();
        const jDayStart = new Date(now.getUTCFullYear(), 0, 0);
        const jDayDiff = now.getDate() - jDayStart.getDate();
        return Math.floor(jDayDiff / _main_js__WEBPACK_IMPORTED_MODULE_0__.MILLISECONDS_TO_DAYS);
    }
    if (typeof mon === 'undefined' ||
        typeof day === 'undefined' ||
        typeof hr === 'undefined' ||
        typeof minute === 'undefined' ||
        typeof sec === 'undefined') {
        throw new Error('Invalid date');
    }
    return (367.0 * year -
        Math.floor(7 * (year + Math.floor((mon + 9) / 12.0)) * 0.25) +
        Math.floor((275 * mon) / 9.0) +
        day +
        1721013.5 +
        ((sec / 60.0 + minute) / 60.0 + hr) / 24.0);
};
/**
 * Calculates the Greenwich Mean Sidereal Time (GMST) for a given date.
 * @param date - The date for which to calculate the GMST.
 * @returns An object containing the GMST value and the Julian date.
 */
function calcGmst(date) {
    const j = jday(date.getUTCFullYear(), date.getUTCMonth() + 1, date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds()) +
        date.getUTCMilliseconds() * _main_js__WEBPACK_IMPORTED_MODULE_0__.MILLISECONDS_TO_DAYS;
    const gmst = _main_js__WEBPACK_IMPORTED_MODULE_0__.Sgp4.gstime(j);
    return { gmst, j };
}
/**
 * Converts ECI coordinates to RAE (Right Ascension, Azimuth, Elevation) coordinates.
 * @variation cached - results are cached
 * @param now - Current date and time.
 * @param eci - ECI coordinates of the satellite.
 * @param sensor - Sensor object containing observer's geodetic coordinates.
 * @returns Object containing azimuth, elevation and range in degrees and kilometers respectively.
 */
function eci2rae(now, eci, sensor) {
    now = new Date(now);
    const { gmst } = calcGmst(now);
    const positionEcf = eci2ecf(eci, gmst);
    const lla = {
        lat: (sensor.lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        lon: (sensor.lon * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD),
        alt: sensor.alt,
    };
    const rae = ecfRad2rae(lla, positionEcf);
    return rae;
}
/**
 * Calculates the inertial azimuth of a satellite given its latitude and inclination.
 * @param lat - The latitude of the satellite in degrees.
 * @param inc - The inclination of the satellite in degrees.
 * @returns The inertial azimuth of the satellite in degrees.
 */
function calcInertAz(lat, inc) {
    const phi = lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    const i = inc * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    const az = Math.asin(Math.cos(i) / Math.cos(phi));
    return (az * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
}
/**
 * Calculates the inclination angle of a satellite from its launch azimuth and latitude.
 * @param lat - The latitude of the observer in degrees.
 * @param az - The launch azimuth angle of the satellite in degrees clockwise from north.
 * @returns The inclination angle of the satellite in degrees.
 */
function calcIncFromAz(lat, az) {
    const phi = lat * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    const beta = az * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    const inc = Math.acos(Math.sin(beta) * Math.cos(phi));
    return (inc * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG);
}
/**
 * Converts Azimuth and Elevation to U and V.
 * Azimuth is the angle off of boresight in the horizontal plane.
 * Elevation is the angle off of boresight in the vertical plane.
 * Cone half angle is the angle of the cone of the radar max field of view.
 * @param az - Azimuth in radians
 * @param el - Elevation in radians
 * @param coneHalfAngle - Cone half angle in radians
 * @returns U and V in radians
 */
function azel2uv(az, el, coneHalfAngle) {
    if (az > coneHalfAngle && az < coneHalfAngle) {
        throw new RangeError(`Azimuth is out of bounds: ${az}`);
    }
    if (el > coneHalfAngle && el < coneHalfAngle) {
        throw new RangeError(`Elevation is out of bounds: ${el}`);
    }
    const alpha = (az / (coneHalfAngle * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG)) * 90;
    const beta = (el / (coneHalfAngle * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG)) * 90;
    const u = Math.sin(alpha);
    let v = -Math.sin(beta);
    v = Object.is(v, -0) ? 0 : v;
    return { u, v };
}
/**
 * Determine azimuth and elevation off of boresight based on sensor orientation and RAE.
 * @param rae Range, Azimuth, Elevation
 * @param sensor Radar sensor object
 * @param face Face number of the sensor
 * @param maxSensorAz Maximum sensor azimuth
 * @returns Azimuth and Elevation off of boresight
 */
function rae2raeOffBoresight(rae, sensor, face, maxSensorAz) {
    let az = (rae.az * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    let el = (rae.el * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    // Correct azimuth for sensor orientation.
    az = az > maxSensorAz * _main_js__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD ? (az - _main_js__WEBPACK_IMPORTED_MODULE_0__.TAU) : az;
    az = (az - sensor.boresightAz[face]);
    el = (el - sensor.boresightEl[face]);
    return { az, el };
}
/**
 * Converts Range Az El to Range U V.
 * @param rae Range, Azimuth, Elevation
 * @param sensor Radar sensor object
 * @param face Face number of the sensor
 * @param maxSensorAz Maximum sensor azimuth
 * @returns Range, U, V
 */
function rae2ruv(rae, sensor, face, maxSensorAz) {
    const { az, el } = rae2raeOffBoresight(rae, sensor, face, maxSensorAz);
    const { u, v } = azel2uv(az, el, sensor.beamwidthRad);
    return { rng: rae.rng, u, v };
}
/**
 * Converts U and V to Azimuth and Elevation off of boresight.
 * @param u The U coordinate.
 * @param v The V coordinate.
 * @param coneHalfAngle The cone half angle of the radar.
 * @returns Azimuth and Elevation off of boresight.
 */
function uv2azel(u, v, coneHalfAngle) {
    if (u > 1 || u < -1) {
        throw new RangeError(`u is out of bounds: ${u}`);
    }
    if (v > 1 || v < -1) {
        throw new RangeError(`v is out of bounds: ${v}`);
    }
    const alpha = Math.asin(u);
    const beta = Math.asin(v);
    const az = ((alpha / 90) * (coneHalfAngle * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG));
    const el = ((beta / 90) * (coneHalfAngle * _main_js__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG));
    return { az, el };
}
/**
 * Converts coordinates from East-North-Up (ENU) to Right-Front-Up (RF) coordinate system.
 * @param enu - The ENU coordinates to be converted.
 * @param enu.x - The east coordinate.
 * @param enu.y - The north coordinate.
 * @param enu.z - The up coordinate.
 * @param az - The azimuth angle in radians.
 * @param el - The elevation angle in radians.
 * @returns The converted RF coordinates.
 */
function enu2rf({ x, y, z }, az, el) {
    const xrf = Math.cos(el) * Math.cos(az) * x - Math.sin(az) * y + Math.sin(el) * Math.cos(az) * z;
    const yrf = Math.cos(el) * Math.sin(az) * x + Math.cos(az) * y + Math.sin(el) * Math.sin(az) * z;
    const zrf = -Math.sin(el) * x + Math.cos(el) * z;
    return {
        x: xrf,
        y: yrf,
        z: zrf,
    };
}


}),
"./src/engine/ootk/src/types/types.ts": 
/*!********************************************!*\
  !*** ./src/engine/ootk/src/types/types.ts ***!
  \********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  PayloadStatus: () => (PayloadStatus),
  SpaceObjectType: () => (SpaceObjectType),
  ZoomValue: () => (ZoomValue)
});
/**
 * @author @thkruz Theodore Kruczek
 * @license AGPL-3.0-or-later
 * @copyright (c) 2025 Kruczek Labs LLC
 *
 * Orbital Object ToolKit is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Orbital Object ToolKit is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * Orbital Object ToolKit. If not, see <http://www.gnu.org/licenses/>.
 */
/**
 * Enum representing different types of objects.
 */
var SpaceObjectType;
(function (SpaceObjectType) {
    SpaceObjectType[SpaceObjectType["UNKNOWN"] = 0] = "UNKNOWN";
    SpaceObjectType[SpaceObjectType["PAYLOAD"] = 1] = "PAYLOAD";
    SpaceObjectType[SpaceObjectType["ROCKET_BODY"] = 2] = "ROCKET_BODY";
    SpaceObjectType[SpaceObjectType["DEBRIS"] = 3] = "DEBRIS";
    SpaceObjectType[SpaceObjectType["SPECIAL"] = 4] = "SPECIAL";
    SpaceObjectType[SpaceObjectType["BALLISTIC_MISSILE"] = 8] = "BALLISTIC_MISSILE";
    SpaceObjectType[SpaceObjectType["STAR"] = 9] = "STAR";
    SpaceObjectType[SpaceObjectType["INTERGOVERNMENTAL_ORGANIZATION"] = 10] = "INTERGOVERNMENTAL_ORGANIZATION";
    SpaceObjectType[SpaceObjectType["SUBORBITAL_PAYLOAD_OPERATOR"] = 11] = "SUBORBITAL_PAYLOAD_OPERATOR";
    SpaceObjectType[SpaceObjectType["PAYLOAD_OWNER"] = 12] = "PAYLOAD_OWNER";
    SpaceObjectType[SpaceObjectType["METEOROLOGICAL_ROCKET_LAUNCH_AGENCY_OR_MANUFACTURER"] = 13] = "METEOROLOGICAL_ROCKET_LAUNCH_AGENCY_OR_MANUFACTURER";
    SpaceObjectType[SpaceObjectType["PAYLOAD_MANUFACTURER"] = 14] = "PAYLOAD_MANUFACTURER";
    SpaceObjectType[SpaceObjectType["LAUNCH_AGENCY"] = 15] = "LAUNCH_AGENCY";
    SpaceObjectType[SpaceObjectType["LAUNCH_SITE"] = 16] = "LAUNCH_SITE";
    SpaceObjectType[SpaceObjectType["LAUNCH_POSITION"] = 17] = "LAUNCH_POSITION";
    SpaceObjectType[SpaceObjectType["LAUNCH_FACILITY"] = 18] = "LAUNCH_FACILITY";
    SpaceObjectType[SpaceObjectType["CONTROL_FACILITY"] = 19] = "CONTROL_FACILITY";
    SpaceObjectType[SpaceObjectType["GROUND_SENSOR_STATION"] = 20] = "GROUND_SENSOR_STATION";
    SpaceObjectType[SpaceObjectType["OPTICAL"] = 21] = "OPTICAL";
    SpaceObjectType[SpaceObjectType["MECHANICAL"] = 22] = "MECHANICAL";
    SpaceObjectType[SpaceObjectType["PHASED_ARRAY_RADAR"] = 23] = "PHASED_ARRAY_RADAR";
    SpaceObjectType[SpaceObjectType["OBSERVER"] = 24] = "OBSERVER";
    SpaceObjectType[SpaceObjectType["BISTATIC_RADIO_TELESCOPE"] = 25] = "BISTATIC_RADIO_TELESCOPE";
    SpaceObjectType[SpaceObjectType["COUNTRY"] = 26] = "COUNTRY";
    SpaceObjectType[SpaceObjectType["LAUNCH_VEHICLE_MANUFACTURER"] = 27] = "LAUNCH_VEHICLE_MANUFACTURER";
    SpaceObjectType[SpaceObjectType["ENGINE_MANUFACTURER"] = 28] = "ENGINE_MANUFACTURER";
    SpaceObjectType[SpaceObjectType["NOTIONAL"] = 29] = "NOTIONAL";
    SpaceObjectType[SpaceObjectType["FRAGMENT"] = 30] = "FRAGMENT";
    SpaceObjectType[SpaceObjectType["SHORT_TERM_FENCE"] = 31] = "SHORT_TERM_FENCE";
    SpaceObjectType[SpaceObjectType["MAX_SPACE_OBJECT_TYPE"] = 32] = "MAX_SPACE_OBJECT_TYPE";
    SpaceObjectType[SpaceObjectType["TERRESTRIAL_PLANET"] = 33] = "TERRESTRIAL_PLANET";
    SpaceObjectType[SpaceObjectType["GAS_GIANT"] = 34] = "GAS_GIANT";
    SpaceObjectType[SpaceObjectType["ICE_GIANT"] = 35] = "ICE_GIANT";
    SpaceObjectType[SpaceObjectType["DWARF_PLANET"] = 36] = "DWARF_PLANET";
    SpaceObjectType[SpaceObjectType["MOON"] = 37] = "MOON";
})(SpaceObjectType || (SpaceObjectType = {}));
var ZoomValue;
(function (ZoomValue) {
    ZoomValue[ZoomValue["LEO"] = 0.45] = "LEO";
    ZoomValue[ZoomValue["GEO"] = 0.82] = "GEO";
    ZoomValue[ZoomValue["MAX"] = 1] = "MAX";
})(ZoomValue || (ZoomValue = {}));
/*
 * + Operational
 * - Nonoperational
 * P Partially Operational
 * Partially fulfilling primary mission or secondary mission(s)
 * B Backup/Standby
 * Previously operational satellite put into reserve status
 * S Spare
 * New satellite awaiting full activation
 * X Extended Mission
 * D Decayed
 * ? Unknown
 */
var PayloadStatus;
(function (PayloadStatus) {
    PayloadStatus["OPERATIONAL"] = "+";
    PayloadStatus["NONOPERATIONAL"] = "-";
    PayloadStatus["PARTIALLY_OPERATIONAL"] = "P";
    PayloadStatus["BACKUP_STANDBY"] = "B";
    PayloadStatus["SPARE"] = "S";
    PayloadStatus["EXTENDED_MISSION"] = "X";
    PayloadStatus["DECAYED"] = "D";
    PayloadStatus["UNKNOWN"] = "?";
})(PayloadStatus || (PayloadStatus = {}));


}),
"./src/engine/ootk/src/utils/constants.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/utils/constants.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DEG2RAD: () => (DEG2RAD),
  MILLISECONDS_PER_DAY: () => (MILLISECONDS_PER_DAY),
  MILLISECONDS_PER_SECOND: () => (MILLISECONDS_PER_SECOND),
  MILLISECONDS_TO_DAYS: () => (MILLISECONDS_TO_DAYS),
  MINUTES_PER_DAY: () => (MINUTES_PER_DAY),
  MS_PER_DAY: () => (MS_PER_DAY),
  PI: () => (PI),
  RAD2DEG: () => (RAD2DEG),
  RADIUS_OF_EARTH: () => (RADIUS_OF_EARTH),
  TAU: () => (TAU),
  angularVelocityOfEarth: () => (angularVelocityOfEarth),
  asec2rad: () => (asec2rad),
  astronomicalUnit: () => (astronomicalUnit),
  cKmPerMs: () => (cKmPerMs),
  cKmPerSec: () => (cKmPerSec),
  cMPerSec: () => (cMPerSec),
  earthGravityParam: () => (earthGravityParam),
  halfPi: () => (halfPi),
  masec2rad: () => (masec2rad),
  msec2sec: () => (msec2sec),
  sec2day: () => (sec2day),
  sec2deg: () => (sec2deg),
  sec2min: () => (sec2min),
  secondsPerDay: () => (secondsPerDay),
  secondsPerSiderealDay: () => (secondsPerSiderealDay),
  secondsPerWeek: () => (secondsPerWeek),
  temp4: () => (temp4),
  ttasec2rad: () => (ttasec2rad),
  x2o3: () => (x2o3)
});
/**
 * Full circle in radians (PI * 2)
 *
 * https://tauday.com/tau-manifesto
 */
const TAU = (2.0 * Math.PI);
/**
 * Represents half of the mathematical constant PI.
 */
const halfPi = (0.5 * Math.PI);
/**
 * Converts degrees to radians.
 */
const DEG2RAD = (Math.PI / 180.0);
/**
 * Converts radians to degrees.
 */
const RAD2DEG = (180.0 / Math.PI);
/**
 * Conversion factor from seconds to degrees.
 */
const sec2deg = (1.0 / 60.0 / 60.0);
/**
 * Conversion factor from seconds to days.
 */
const sec2day = sec2deg / 24.0;
/**
 * Conversion factor from arcseconds to radians.
 */
const asec2rad = (sec2deg * DEG2RAD);
/**
 * Convert ten-thousandths of an arcsecond to radians.
 */
const ttasec2rad = (asec2rad / 10000.0);
/**
 * Convert milliarcseconds to radians.
 */
const masec2rad = (asec2rad / 1000.0);
/**
 * The angular velocity of the Earth in radians per second.
 */
const angularVelocityOfEarth = 7.292115e-5;
/**
 * Astronomical unit in kilometers.
 */
const astronomicalUnit = 149597870.0;
// / Convert milliseconds to seconds.
const msec2sec = 1e-3;
// / Speed of light.
const cMPerSec = 299792458;
const cKmPerSec = 299792458 / 1000;
const cKmPerMs = 299792458 / 1000 / 1000;
// / Milliseconds per day.
const MS_PER_DAY = 86400000;
// / Seconds per day.
const secondsPerDay = 86400.0;
// / Convert seconds to minutes.
const sec2min = (1.0 / 60.0);
// / Seconds per sidereal day.
const secondsPerSiderealDay = 86164.0905;
// / Seconds per week.
const secondsPerWeek = secondsPerDay * 7.0;
/**
 * Half the number of radians in a circle.
 */
const PI = Math.PI;
const x2o3 = 2.0 / 3.0;
const temp4 = 1.5e-12;
/**
 * The number of minutes in a day.
 */
const MINUTES_PER_DAY = 1440;
/**
 * The number of milliseconds in a day.
 */
const MILLISECONDS_TO_DAYS = 1.15741e-8;
/**
 * The number of milliseconds in a day.
 */
const MILLISECONDS_PER_DAY = 1000 * 60 * 60 * 24;
/**
 * The number of milliseconds in a second.
 */
const MILLISECONDS_PER_SECOND = 1000;
const RADIUS_OF_EARTH = 6371; // Radius of Earth in kilometers
const earthGravityParam = 398600.4415;


}),
"./src/engine/ootk/src/utils/create-covariance-from-tle.ts": 
/*!*****************************************************************!*\
  !*** ./src/engine/ootk/src/utils/create-covariance-from-tle.ts ***!
  \*****************************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  createCovarianceFromTle: () => (createCovarianceFromTle),
  createSampleCovarianceFromTle: () => (createSampleCovarianceFromTle)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");

/**
 * Creates a 6x6 state covariance matrix from a TLE
 * @param tleLine1 The first line of the TLE
 * @param tleLine2 The second line of the TLE
 * @param frame The covariance frame (CovarianceFrame.ECI or CovarianceFrame.RIC)
 * @param sigmaScale Scaling factor for the sigmas (default: 1.0)
 * @returns A StateCovariance object containing the 6x6 covariance matrix
 */
function createCovarianceFromTle(tleLine1, tleLine2, frame = _main_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceFrame.RIC, sigmaScale = 1.0) {
    // Parse the TLE and get the state vector
    const tle = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Tle(tleLine1, tleLine2);
    const temeState = tle.state;
    const j2000State = temeState.toJ2000();
    // eslint-disable-next-line no-console
    console.log('J2000 State:', j2000State);
    /*
     * Determine appropriate sigma values based on TLE
     * These are rough estimates and can be adjusted based on your needs
     */
    const positionSigma = 1.0 * sigmaScale; // km
    const velocitySigma = 0.001 * sigmaScale; // km/s
    // Create a diagonal covariance matrix with sigma values
    const sigmas = [
        positionSigma, positionSigma, positionSigma,
        velocitySigma, velocitySigma, velocitySigma,
    ];
    // Generate the covariance matrix in the desired frame
    return _main_js__WEBPACK_IMPORTED_MODULE_0__.StateCovariance.fromSigmas(sigmas, frame);
}
/**
 * Creates a sample-based covariance from a TLE with more realistic uncertainties
 * @param tleLine1 The first line of the TLE
 * @param tleLine2 The second line of the TLE
 * @param frame The covariance frame (CovarianceFrame.ECI or CovarianceFrame.RIC)
 * @returns A StateCovariance object
 */
function createSampleCovarianceFromTle(tleLine1, tleLine2, frame = _main_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceFrame.RIC) {
    // Parse the TLE and get the state vector
    const tle = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Tle(tleLine1, tleLine2);
    const temeState = tle.state;
    const j2000State = temeState.toJ2000();
    /*
     * Create initial covariance with basic sigma values
     * Position uncertainties are higher in-track than radial/cross-track
     * for most space catalog objects
     */
    const sigmas = frame === _main_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceFrame.RIC
        ? [0.12, 1.0, 0.1, 0.00012, 0.001, 0.0001] // RIC frame: [R,I,C,Rdot,Idot,Cdot]
        : [0.6, 0.6, 0.6, 0.0006, 0.0006, 0.0006]; // ECI frame: [x,y,z,vx,vy,vz]
    const covariance = _main_js__WEBPACK_IMPORTED_MODULE_0__.StateCovariance.fromSigmas(sigmas, frame);
    // Create a covariance sample that will be used to generate a more realistic covariance
    const sample = new _main_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceSample(j2000State, covariance, tle);
    // Return the desample in the appropriate frame
    return frame === _main_js__WEBPACK_IMPORTED_MODULE_0__.CovarianceFrame.RIC ? sample.desampleRIC() : sample.desampleJ2000();
}


}),
"./src/engine/ootk/src/utils/functions.ts": 
/*!************************************************!*\
  !*** ./src/engine/ootk/src/utils/functions.ts ***!
  \************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  acoth: () => (acoth),
  acsch: () => (acsch),
  angularDiameter: () => (angularDiameter),
  angularDistance: () => (angularDistance),
  array2d: () => (array2d),
  asech: () => (asech),
  clamp: () => (clamp),
  concat: () => (concat),
  copySign: () => (copySign),
  covariance: () => (covariance),
  createVec: () => (createVec),
  csch: () => (csch),
  derivative: () => (derivative),
  dopplerFactor: () => (dopplerFactor),
  evalPoly: () => (evalPoly),
  factorial: () => (factorial),
  gamma: () => (gamma),
  getDayOfYear: () => (getDayOfYear),
  isLeapYear: () => (isLeapYear),
  linearInterpolate: () => (linearInterpolate),
  log10: () => (log10),
  matchHalfPlane: () => (matchHalfPlane),
  mean: () => (mean),
  newtonM: () => (newtonM),
  newtonNu: () => (newtonNu),
  sech: () => (sech),
  sign: () => (sign),
  spaceObjType2Str: () => (spaceObjType2Str),
  std: () => (std),
  toPrecision: () => (toPrecision),
  wrapAngle: () => (wrapAngle)
});
/* ESM import */var _enums_AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../enums/AngularDiameterMethod.js */ "./src/engine/ootk/src/enums/AngularDiameterMethod.ts");
/* ESM import */var _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../enums/AngularDistanceMethod.js */ "./src/engine/ootk/src/enums/AngularDistanceMethod.ts");
/* ESM import */var _types_types_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../types/types.js */ "./src/engine/ootk/src/types/types.ts");
/* ESM import */var _constants_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* eslint-disable require-jsdoc */




/**
 * Calculates the factorial of a given number.
 * @param n - The number to calculate the factorial for.
 * @returns The factorial of the given number.
 */
function factorial(n) {
    const nAbs = Math.abs(n);
    let result = 1;
    for (let i = 2; i <= nAbs; i++) {
        result *= i;
    }
    return result;
}
/**
 * Calculates the base 10 logarithm of a number.
 * @param x - The number to calculate the logarithm for.
 * @returns The base 10 logarithm of the input number.
 */
function log10(x) {
    return Math.log(x) / Math.LN10;
}
/**
 * Calculates the hyperbolic secant of a number.
 * @param x - The number to calculate the hyperbolic secant of.
 * @returns The hyperbolic secant of the given number.
 */
function sech(x) {
    return 1 / Math.cosh(x);
}
/**
 * Calculates the hyperbolic cosecant of a number.
 * @param x - The number for which to calculate the hyperbolic cosecant.
 * @returns The hyperbolic cosecant of the given number.
 */
function csch(x) {
    return 1 / Math.sinh(x);
}
/**
 * Returns the inverse hyperbolic cosecant of a number.
 * @param x - The number to calculate the inverse hyperbolic cosecant of.
 * @returns The inverse hyperbolic cosecant of the given number.
 */
function acsch(x) {
    return Math.log(1 / x + Math.sqrt(1 / (x * x) + 1));
}
/**
 * Calculates the inverse hyperbolic secant (asech) of a number.
 * @param x - The number to calculate the inverse hyperbolic secant of.
 * @returns The inverse hyperbolic secant of the given number.
 */
function asech(x) {
    return Math.log(1 / x + Math.sqrt(1 / (x * x) - 1));
}
/**
 * Calculates the inverse hyperbolic cotangent (acoth) of a number.
 * @param x - The number to calculate the acoth of.
 * @returns The inverse hyperbolic cotangent of the given number.
 */
function acoth(x) {
    return 0.5 * Math.log((x + 1) / (x - 1));
}
/**
 * Copies the sign of the second number to the first number.
 * @param mag - The magnitude of the number.
 * @param sgn - The sign of the number.
 * @returns The number with the magnitude of `mag` and the sign of `sgn`.
 */
function copySign(mag, sgn) {
    return Math.abs(mag) * Math.sign(sgn);
}
/**
 * Evaluates a polynomial function at a given value.
 * @param x - The value at which to evaluate the polynomial.
 * @param coeffs - The coefficients of the polynomial.
 * @returns The result of evaluating the polynomial at the given value.
 */
function evalPoly(x, coeffs) {
    let result = coeffs[0];
    for (let i = 1; i < coeffs.length; i++) {
        result = result * x + coeffs[i];
    }
    return result;
}
/**
 * Concatenates two Float64Arrays into a new Float64Array.
 * @param a - The first Float64Array.
 * @param b - The second Float64Array.
 * @returns A new Float64Array containing the concatenated values of `a` and `b`.
 */
function concat(a, b) {
    const result = new Float64Array(a.length + b.length);
    result.set(a);
    result.set(b, a.length);
    return result;
}
/**
 * Calculates the angle in the half-plane that best matches the given angle.
 * @param angle - The angle to be matched.
 * @param match - The angle to be matched against.
 * @returns The angle in the half-plane that best matches the given angle.
 */
function matchHalfPlane(angle, match) {
    const a1 = angle;
    const a2 = 2 * Math.PI - angle;
    const d1 = Math.atan2(Math.sin(a1 - match), Math.cos(a1 - match));
    const d2 = Math.atan2(Math.sin(a2 - match), Math.cos(a2 - match));
    return Math.abs(d1) < Math.abs(d2) ? a1 : a2;
}
/**
 * Wraps an angle to the range [-π, π].
 * @param theta - The angle to wrap.
 * @returns The wrapped angle.
 */
function wrapAngle(theta) {
    const result = ((theta + Math.PI) % (2 * Math.PI)) - Math.PI;
    if (result === -Math.PI) {
        return Math.PI;
    }
    return result;
}
/**
 * Calculates the angular distance between two points on a sphere using the cosine formula.
 * @param lam1 - The longitude of the first point in radians.
 * @param phi1 - The latitude of the first point in radians.
 * @param lam2 - The longitude of the second point in radians.
 * @param phi2 - The latitude of the second point in radians.
 * @returns The angular distance between the two points in radians.
 */
function angularDistanceCosine_(lam1, phi1, lam2, phi2) {
    const a = Math.sin(phi1) * Math.sin(phi2);
    const b = Math.cos(phi1) * Math.cos(phi2) * Math.cos(lam2 - lam1);
    return Math.acos(a + b);
}
/**
 * Calculates the angular distance between two points on a sphere using the Haversine formula.
 * @param lam1 - The longitude of the first point in radians.
 * @param phi1 - The latitude of the first point in radians.
 * @param lam2 - The longitude of the second point in radians.
 * @param phi2 - The latitude of the second point in radians.
 * @returns The angular distance between the two points in radians.
 */
function angularDistanceHaversine_(lam1, phi1, lam2, phi2) {
    const dlam = lam2 - lam1;
    const dphi = phi2 - phi1;
    const sdlam = Math.sin(0.5 * dlam);
    const sdphi = Math.sin(0.5 * dphi);
    const a = sdphi * sdphi + Math.cos(phi1) * Math.cos(phi2) * sdlam * sdlam;
    return 2.0 * Math.asin(Math.min(1.0, Math.sqrt(a)));
}
/**
 * Calculates the angular distance between two points on a sphere.
 * @param lam1 The longitude of the first point.
 * @param phi1 The latitude of the first point.
 * @param lam2 The longitude of the second point.
 * @param phi2 The latitude of the second point.
 * @param method The method to use for calculating the angular distance. Defaults to AngularDistanceMethod.Cosine.
 * @returns The angular distance between the two points.
 * @throws Error if an invalid angular distance method is provided.
 */
function angularDistance(lam1, phi1, lam2, phi2, method = _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Cosine) {
    switch (method) {
        case _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Cosine:
            return angularDistanceCosine_(lam1, phi1, lam2, phi2);
        case _enums_AngularDistanceMethod_js__WEBPACK_IMPORTED_MODULE_1__.AngularDistanceMethod.Haversine:
            return angularDistanceHaversine_(lam1, phi1, lam2, phi2);
        default:
            throw new Error('Invalid angular distance method.');
    }
}
/**
 * Calculates the angular diameter of an object.
 * @param diameter - The diameter of the object.
 * @param distance - The distance to the object.
 * @param method - The method used to calculate the angular diameter. Defaults to AngularDiameterMethod.Sphere.
 * @returns The angular diameter of the object.
 * @throws Error if an invalid angular diameter method is provided.
 */
function angularDiameter(diameter, distance, method = _enums_AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Sphere) {
    switch (method) {
        case _enums_AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Circle:
            return 2 * Math.atan(diameter / (2 * distance));
        case _enums_AngularDiameterMethod_js__WEBPACK_IMPORTED_MODULE_0__.AngularDiameterMethod.Sphere:
            return 2 * Math.asin(diameter / (2 * distance));
        default:
            throw new Error('Invalid angular diameter method.');
    }
}
/**
 * Performs linear interpolation between two points.
 * @param x - The x-coordinate to interpolate.
 * @param x0 - The x-coordinate of the first point.
 * @param y0 - The y-coordinate of the first point.
 * @param x1 - The x-coordinate of the second point.
 * @param y1 - The y-coordinate of the second point.
 * @returns The interpolated y-coordinate corresponding to the given x-coordinate.
 */
function linearInterpolate(x, x0, y0, x1, y1) {
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}
/**
 * Calculates the mean value of an array of numbers.
 * @param values - The array of numbers.
 * @returns The mean value of the numbers.
 */
function mean(values) {
    const n = values.length;
    let sum = 0.0;
    for (const v of values) {
        sum += v;
    }
    return sum / n;
}
/**
 * Calculates the standard deviation of an array of numbers.
 * @param values - The array of numbers.
 * @param isSample - Optional. Specifies whether the array represents a sample. Default is false.
 * @returns The standard deviation of the array.
 */
function std(values, isSample = false) {
    const mu = mean(values);
    const n = values.length;
    let sum = 0.0;
    for (const v of values) {
        const sub = v - mu;
        sum += sub * sub;
    }
    const m = isSample ? 1 : 0;
    return Math.sqrt((1.0 / (n - m)) * sum);
}
/**
 * Calculates the covariance between two arrays.
 * @param a - The first array.
 * @param b - The second array.
 * @param isSample - Optional. Specifies whether the arrays represent a sample. Default is false.
 * @returns The covariance between the two arrays.
 */
function covariance(a, b, isSample = false) {
    const n = a.length;
    const am = mean(a);
    const bm = mean(b);
    let result = 0.0;
    for (let i = 0; i < n; i++) {
        result += (a[i] - am) * (b[i] - bm);
    }
    const m = isSample ? 1 : 0;
    return result / (n - m);
}
/**
 * Calculates the gamma function of a number.
 * @param n - The input number.
 * @returns The gamma function value.
 */
function gamma(n) {
    return factorial(n - 1);
}
/**
 * Calculates the eccentric anomaly (e0) and true anomaly (nu) using Newton's method
 * for a given eccentricity (ecc) and mean anomaly (m).
 * @param ecc - The eccentricity of the orbit.
 * @param m - The mean anomaly.
 * @returns An object containing the eccentric anomaly (e0) and true anomaly (nu).
 */
function newtonM(ecc, m) {
    const numiter = 50;
    const small = 1e-8;
    let e0;
    let nu;
    if (ecc > small) {
        if ((m < 0.0 && m > -Math.PI) || m > Math.PI) {
            e0 = m - ecc;
        }
        else {
            e0 = m + ecc;
        }
        let ktr = 1;
        let e1 = e0 + (m - e0 + ecc * Math.sin(e0)) / (1.0 - ecc * Math.cos(e0));
        while (Math.abs(e1 - e0) > small && ktr <= numiter) {
            ktr++;
            e0 = e1;
            e1 = e0 + (m - e0 + ecc * Math.sin(e0)) / (1.0 - ecc * Math.cos(e0));
        }
        const sinv = (Math.sqrt(1.0 - ecc * ecc) * Math.sin(e1)) / (1.0 - ecc * Math.cos(e1));
        const cosv = (Math.cos(e1) - ecc) / (1.0 - ecc * Math.cos(e1));
        nu = Math.atan2(sinv, cosv);
    }
    else {
        nu = m;
        e0 = m;
    }
    return { e0, nu };
}
/**
 * Calculates the eccentric anomaly (e0) and mean anomaly (m) using Newton's method
 * for a given eccentricity (ecc) and true anomaly (nu).
 * @param ecc - The eccentricity of the orbit.
 * @param nu - The true anomaly.
 * @returns An object containing the calculated eccentric anomaly (e0) and mean anomaly (m).
 */
function newtonNu(ecc, nu) {
    const small = 1e-8;
    let e0 = 0.0;
    let m = 0.0;
    if (Math.abs(ecc) < small) {
        m = nu;
        e0 = nu;
    }
    else if (ecc < 1.0 - small) {
        const sine = (Math.sqrt(1.0 - ecc * ecc) * Math.sin(nu)) / (1.0 + ecc * Math.cos(nu));
        const cose = (ecc + Math.cos(nu)) / (1.0 + ecc * Math.cos(nu));
        e0 = Math.atan2(sine, cose);
        m = e0 - ecc * Math.sin(e0);
    }
    if (ecc < 1.0) {
        m -= Math.floor(m / (2 * Math.PI)) * (2 * Math.PI);
        if (m < 0.0) {
            m += 2.0 * Math.PI;
        }
        e0 -= Math.floor(e0 / (2 * Math.PI)) * (2 * Math.PI);
    }
    return { e0, m: m };
}
/**
 * Creates a 2D array with the specified number of rows and columns, filled with the same given value.
 * @template T The type of elements in the array.
 * @param rows The number of rows in the 2D array.
 * @param columns The number of columns in the 2D array.
 * @param value The value to fill the array with.
 * @returns The 2D array with the specified number of rows and columns, filled with the given value.
 */
function array2d(rows, columns, value) {
    const output = [];
    for (let i = 0; i < rows; i++) {
        output.push(Array(columns).fill(value));
    }
    return output;
}
/**
 * Clamps a number between a minimum and maximum value.
 * @param x The number to clamp.
 * @param min The minimum value.
 * @param max The maximum value.
 * @returns The clamped number.
 */
function clamp(x, min, max) {
    return Math.max(min, Math.min(x, max));
}
/**
 * Determines whether a given year is a leap year.
 * @param dateIn The date to check.
 * @returns `true` if the year is a leap year, `false` otherwise.
 */
function isLeapYear(dateIn) {
    const year = dateIn.getUTCFullYear();
    if ((year & 3) !== 0) {
        return false;
    }
    return year % 100 !== 0 || year % 400 === 0;
}
/**
 * Calculates the day of the year for a given date.
 * If no date is provided, the current date is used.
 *
 * This is sometimes referred to as the Jday, but is
 * very different from the Julian day used in astronomy.
 * @param date - The date for which to calculate the day of the year.
 * @returns The day of the year as a number.
 */
function getDayOfYear(date = new Date()) {
    const dayCount = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    const mn = date.getUTCMonth();
    const dn = date.getUTCDate();
    let dayOfYear = dayCount[mn] + dn;
    if (mn > 1 && isLeapYear(date)) {
        dayOfYear++;
    }
    return dayOfYear;
}
/**
 * Rounds a number to a specified number of decimal places.
 * @param value - The number to round.
 * @param places - The number of decimal places to round to.
 * @returns The rounded number.
 */
function toPrecision(value, places) {
    return Number(value.toFixed(places));
}
/**
 * Returns the sign of a number.
 * @param value - The number to determine the sign of.
 * @returns 1 if the number is positive, -1 if the number is negative.
 */
function sign(value) {
    return value >= 0 ? 1 : -1;
}
const spaceObjTypeStrMap_ = {
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.UNKNOWN]: 'Unknown',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.PAYLOAD]: 'Payload',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.ROCKET_BODY]: 'Rocket Body',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.DEBRIS]: 'Debris',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.SPECIAL]: 'Special',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.BALLISTIC_MISSILE]: 'Ballistic Missile',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.STAR]: 'Star',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.INTERGOVERNMENTAL_ORGANIZATION]: 'Intergovernmental Organization',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.SUBORBITAL_PAYLOAD_OPERATOR]: 'Suborbital Payload Operator',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.PAYLOAD_OWNER]: 'Payload Owner',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.METEOROLOGICAL_ROCKET_LAUNCH_AGENCY_OR_MANUFACTURER]: 'Meteorological Rocket Launch Agency or Manufacturer',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.PAYLOAD_MANUFACTURER]: 'Payload Manufacturer',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.LAUNCH_AGENCY]: 'Launch Agency',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.LAUNCH_SITE]: 'Launch Site',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.LAUNCH_POSITION]: 'Launch Position',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.LAUNCH_FACILITY]: 'Launch Facility',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.CONTROL_FACILITY]: 'Control Facility',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.GROUND_SENSOR_STATION]: 'Ground Sensor Station',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.OPTICAL]: 'Optical',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.MECHANICAL]: 'Mechanical',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.PHASED_ARRAY_RADAR]: 'Phased Array Radar',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.OBSERVER]: 'Observer',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.BISTATIC_RADIO_TELESCOPE]: 'Bi-static Radio Telescope',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.COUNTRY]: 'Country',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.LAUNCH_VEHICLE_MANUFACTURER]: 'Launch Vehicle Manufacturer',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.ENGINE_MANUFACTURER]: 'Engine Manufacturer',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.NOTIONAL]: 'Notional',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.FRAGMENT]: 'Fragment',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.SHORT_TERM_FENCE]: 'Short Term Fence',
    [_types_types_js__WEBPACK_IMPORTED_MODULE_2__.SpaceObjectType.MAX_SPACE_OBJECT_TYPE]: 'Max Space Object Type',
};
/**
 * Converts a SpaceObjectType to a string representation.
 * @param spaceObjType - The SpaceObjectType to convert.
 * @returns The string representation of the SpaceObjectType.
 */
const spaceObjType2Str = (spaceObjType) => spaceObjTypeStrMap_[spaceObjType] || 'Unknown';
/**
 * Calculates the Doppler factor for a given location, position, and velocity.
 * The Doppler factor is a measure of the change in frequency or wavelength of a wave
 * as observed by an observer moving relative to the source of the wave.
 * @param location - The location vector of the observer.
 * @param position - The position vector of the source.
 * @param velocity - The velocity vector of the source.
 * @returns The calculated Doppler factor.
 */
const dopplerFactor = (location, position, velocity) => {
    const range = {
        x: position.x - location.x,
        y: position.y - location.y,
        z: position.z - location.z,
    };
    const distance = Math.sqrt(range.x ** 2 + range.y ** 2 + range.z ** 2);
    const rangeVel = {
        x: velocity.x + _constants_js__WEBPACK_IMPORTED_MODULE_3__.angularVelocityOfEarth * location.y,
        y: velocity.y - _constants_js__WEBPACK_IMPORTED_MODULE_3__.angularVelocityOfEarth * location.x,
        z: velocity.z,
    };
    const rangeRate = (range.x * rangeVel.x + range.y * rangeVel.y + range.z * rangeVel.z) / distance;
    const dopplerFactor = 1 - rangeRate / _constants_js__WEBPACK_IMPORTED_MODULE_3__.cKmPerSec;
    return dopplerFactor;
};
/**
 * Creates an array of numbers from start to stop (inclusive) with the specified step.
 * @param start The starting number.
 * @param stop The ending number.
 * @param step The step value.
 * @returns An array of numbers.
 */
function createVec(start, stop, step) {
    const array = [];
    for (let i = start; i <= stop; i += step) {
        array.push(i);
    }
    return array;
}
/**
 * Calculates the derivative of a differentiable function.
 * @param f The differentiable function.
 * @param h The step size for numerical differentiation. Default value is 1e-3.
 * @returns The derivative function.
 */
function derivative(f, h = 1e-3) {
    /**
     * @param x The value at which to calculate the derivative.
     * @returns The derivative of the function at the given value.
     */
    function df(x) {
        const hh = h * 0.5;
        return (f(x + hh) - f(x - hh)) / h;
    }
    return df;
}


}),
"./src/engine/ootk/src/utils/index.ts": 
/*!********************************************!*\
  !*** ./src/engine/ootk/src/utils/index.ts ***!
  \********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DEG2RAD: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.DEG2RAD),
  MILLISECONDS_PER_DAY: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.MILLISECONDS_PER_DAY),
  MILLISECONDS_PER_SECOND: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.MILLISECONDS_PER_SECOND),
  MILLISECONDS_TO_DAYS: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.MILLISECONDS_TO_DAYS),
  MINUTES_PER_DAY: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.MINUTES_PER_DAY),
  MS_PER_DAY: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.MS_PER_DAY),
  PI: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.PI),
  RAD2DEG: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.RAD2DEG),
  RADIUS_OF_EARTH: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.RADIUS_OF_EARTH),
  TAU: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.TAU),
  acoth: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.acoth),
  acsch: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.acsch),
  angularDiameter: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.angularDiameter),
  angularDistance: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.angularDistance),
  angularVelocityOfEarth: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.angularVelocityOfEarth),
  array2d: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.array2d),
  asec2rad: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.asec2rad),
  asech: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.asech),
  astronomicalUnit: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.astronomicalUnit),
  cKmPerMs: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.cKmPerMs),
  cKmPerSec: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.cKmPerSec),
  cMPerSec: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.cMPerSec),
  clamp: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.clamp),
  concat: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.concat),
  copySign: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.copySign),
  covariance: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.covariance),
  createCovarianceFromTle: () => (/* reexport safe */ _create_covariance_from_tle_js__WEBPACK_IMPORTED_MODULE_4__.createCovarianceFromTle),
  createSampleCovarianceFromTle: () => (/* reexport safe */ _create_covariance_from_tle_js__WEBPACK_IMPORTED_MODULE_4__.createSampleCovarianceFromTle),
  createVec: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.createVec),
  csch: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.csch),
  derivative: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.derivative),
  dopplerFactor: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.dopplerFactor),
  earthGravityParam: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.earthGravityParam),
  evalPoly: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.evalPoly),
  factorial: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.factorial),
  gamma: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.gamma),
  getDayOfYear: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.getDayOfYear),
  halfPi: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.halfPi),
  isLeapYear: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.isLeapYear),
  jacobian: () => (/* reexport safe */ _jacobian_js__WEBPACK_IMPORTED_MODULE_3__.jacobian),
  linearDistance: () => (/* reexport safe */ _linearDistance_js__WEBPACK_IMPORTED_MODULE_2__.linearDistance),
  linearInterpolate: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.linearInterpolate),
  log10: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.log10),
  masec2rad: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.masec2rad),
  matchHalfPlane: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.matchHalfPlane),
  mean: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.mean),
  msec2sec: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.msec2sec),
  newtonM: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.newtonM),
  newtonNu: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.newtonNu),
  sec2day: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.sec2day),
  sec2deg: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.sec2deg),
  sec2min: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.sec2min),
  sech: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.sech),
  secondsPerDay: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.secondsPerDay),
  secondsPerSiderealDay: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.secondsPerSiderealDay),
  secondsPerWeek: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.secondsPerWeek),
  sign: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.sign),
  spaceObjType2Str: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.spaceObjType2Str),
  std: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.std),
  temp4: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.temp4),
  toPrecision: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.toPrecision),
  ttasec2rad: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.ttasec2rad),
  wrapAngle: () => (/* reexport safe */ _functions_js__WEBPACK_IMPORTED_MODULE_0__.wrapAngle),
  x2o3: () => (/* reexport safe */ _constants_js__WEBPACK_IMPORTED_MODULE_1__.x2o3)
});
/* ESM import */var _functions_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ./functions.js */ "./src/engine/ootk/src/utils/functions.ts");
/* ESM import */var _constants_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./constants.js */ "./src/engine/ootk/src/utils/constants.ts");
/* ESM import */var _linearDistance_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ./linearDistance.js */ "./src/engine/ootk/src/utils/linearDistance.ts");
/* ESM import */var _jacobian_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./jacobian.js */ "./src/engine/ootk/src/utils/jacobian.ts");
/* ESM import */var _create_covariance_from_tle_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./create-covariance-from-tle.js */ "./src/engine/ootk/src/utils/create-covariance-from-tle.ts");







}),
"./src/engine/ootk/src/utils/jacobian.ts": 
/*!***********************************************!*\
  !*** ./src/engine/ootk/src/utils/jacobian.ts ***!
  \***********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  jacobian: () => (jacobian)
});
/* ESM import */var _main_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! ../main.js */ "./src/engine/ootk/src/main.ts");

/**
 * Calculates the Jacobian matrix of a given Jacobian function.
 *
 * The function calculates how small perturbations in each input variable affect all output variables,
 * using a second-order accurate central difference approximation.
 *
 * In orbital mechanics applications, this matrix is essential for solving complex problems like
 * orbit transfers, trajectory optimization, and precise orbital determination.
 * @param f The Jacobian function.
 * @param m The number of rows in the Jacobian matrix.
 * @param x0 The initial values of the variables.
 * @param step The step size for numerical differentiation (default: 1e-5).
 * @returns The Jacobian matrix.
 */
const jacobian = (f, m, x0, step = 0.00001) => {
    const n = x0.length;
    const j = (0,_main_js__WEBPACK_IMPORTED_MODULE_0__.array2d)(m, n, 0);
    const h = 0.5 * step;
    for (let k = 0; k < n; k++) {
        const xp = x0.slice();
        xp[k] += h;
        const fp = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(f(xp));
        const xm = x0.slice();
        xm[k] -= h;
        const fm = new _main_js__WEBPACK_IMPORTED_MODULE_0__.Vector(f(xm));
        const cd = fp.subtract(fm).scale(1 / step);
        for (let i = 0; i < m; i++) {
            j[i][k] = cd.elements[i];
        }
    }
    return new _main_js__WEBPACK_IMPORTED_MODULE_0__.Matrix(j);
};


}),
"./src/engine/ootk/src/utils/linearDistance.ts": 
/*!*****************************************************!*\
  !*** ./src/engine/ootk/src/utils/linearDistance.ts ***!
  \*****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  linearDistance: () => (linearDistance)
});
/**
 * Calculates the linear distance between two points in three-dimensional space.
 * @param pos1 The first position.
 * @param pos2 The second position.
 * @returns The linear distance between the two positions in kilometers.
 */
function linearDistance(pos1, pos2) {
    return Math.sqrt((pos1.x - pos2.x) ** 2 + (pos1.y - pos2.y) ** 2 + (pos1.z - pos2.z) ** 2);
}


}),
"./src/engine/utils/constants.ts": 
/*!***************************************!*\
  !*** ./src/engine/utils/constants.ts ***!
  \***************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  DISTANCE_TO_SUN: () => (DISTANCE_TO_SUN),
  EARTHS_GRAV_CONST: () => (EARTHS_GRAV_CONST),
  EARTH_OBLIQUITY_DEGREES: () => (EARTH_OBLIQUITY_DEGREES),
  EARTH_OBLIQUITY_RADIANS: () => (EARTH_OBLIQUITY_RADIANS),
  GROUND_BUFFER_DISTANCE: () => (GROUND_BUFFER_DISTANCE),
  MASS_OF_EARTH: () => (MASS_OF_EARTH),
  MOON_SCALAR_DISTANCE: () => (MOON_SCALAR_DISTANCE),
  PLANETARIUM_DIST: () => (PLANETARIUM_DIST),
  RADIUS_OF_DRAW_MOON: () => (RADIUS_OF_DRAW_MOON),
  RADIUS_OF_DRAW_SUN: () => (RADIUS_OF_DRAW_SUN),
  RADIUS_OF_EARTH: () => (RADIUS_OF_EARTH),
  RADIUS_OF_SUN: () => (RADIUS_OF_SUN),
  STAR_DISTANCE: () => (STAR_DISTANCE),
  SUN_SCALAR_DISTANCE: () => (SUN_SCALAR_DISTANCE),
  ZOOM_EXP: () => (ZOOM_EXP)
});
/**
 * /////////////////////////////////////////////////////////////////////////////
 *
 * The file constants.ts contains a set of constants used in the KeepTrack application.
 * It defines constants related to zooming, radians, degrees, milliseconds, Earth's gravitational
 * constant, minutes, distances of the Sun and Moon from the Earth, and more.
 * https://github.com/Pownkumar1234/powney---3d-satelite-toolkit
 *
 * @Copyright (C) 2025 Kruczek Labs LLC
 *
 * KeepTrack is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * KeepTrack is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with
 * KeepTrack. If not, see <http://www.gnu.org/licenses/>.
 *
 * /////////////////////////////////////////////////////////////////////////////
 */
/**
 * Exponent used for calculating zoom distance.
 */
const ZOOM_EXP = 3.8;
/**
 * Earths gravitational constant.
 */
const EARTHS_GRAV_CONST = 6.6725985e-11;
/**
 * The mass of the Earth in kilograms.
 */
const MASS_OF_EARTH = 5.97378250603408e24;
/**
 * The offset away from the earth when drawing the camera in Planetarium or Astronomy modes.
 */
const PLANETARIUM_DIST = 3;
/**
 * The radius used to draw the Sun in kilometers.
 * This value is used to scale the Sun's size when rendering it.
 */
const RADIUS_OF_DRAW_SUN = 9000;
/**
 * The scalar distance used to draw the Sun.
 * This value is used to scale the Sun's size and distance from the Earth.
 */
const SUN_SCALAR_DISTANCE = 250000;
/**
 * The radius used to draw the Moon in kilometers.
 * This value is used to scale the Moon's size when rendering it.
 */
const RADIUS_OF_DRAW_MOON = 4000;
/**
 * The scalar distance used to draw the Moon.
 * This value is used to scale the Moon's size and distance from the Earth.
 */
const MOON_SCALAR_DISTANCE = 200000;
/**
 * Radius of the Earth in kilometers.
 */
const RADIUS_OF_EARTH = 6371; // Radius of Earth in kilometers
/**
 * Distance objects are placed above earth to avoid z-buffer fighting
 */
const GROUND_BUFFER_DISTANCE = 2.5;
/**
 * Radius of the Sun in kilometers
 */
const RADIUS_OF_SUN = 695700;
/**
 * Artificial Star Distance - Lower number Reduces webgl depth buffer
 */
const STAR_DISTANCE = 250000;
/**
 * Distance from Earth to the Sun in kilometers
 */
const DISTANCE_TO_SUN = 149597870; // Distance from Earth to the Sun in kilometers
/**
 * Earth's Obliquity in degrees
 */
const EARTH_OBLIQUITY_DEGREES = 23.438480461241912;
/**
 * Earth's Obliquity in radians
 */
const EARTH_OBLIQUITY_RADIANS = EARTH_OBLIQUITY_DEGREES * (Math.PI / 180);


}),
"./src/engine/utils/external/meuusjs.ts": 
/*!**********************************************!*\
  !*** ./src/engine/utils/external/meuusjs.ts ***!
  \**********************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  A: () => (A)
});
/* eslint-disable */
// @ts-nocheck
const A = {
    JMod: 2400000.5,
    J2000: 2451545,
    J1900: 2415020,
    B1900: 2415020.3135,
    B1950: 2433282.4235,
    JulianYear: 365.25,
    JulianCentury: 36525,
    BesselianYear: 365.2421988,
    AU: 149597870,
};
A.EclCoord = function (a, b, c) {
    if (isNaN(a) || isNaN(b))
        throw Error('Invalid EclCoord object: (' + a + ', ' + b + ')');
    this.lat = a;
    this.lng = b;
    void 0 !== c && (this.h = c);
};
A.EclCoord.prototype = {
    toWgs84String: function () {
        return A.Math.formatNum((180 * this.lat) / Math.PI) + ', ' + A.Math.formatNum((180 * -this.lng) / Math.PI);
    },
};
A.EclCoordfromWgs84 = function (a, b, c) {
    return new A.EclCoord((a * Math.PI) / 180, (-b * Math.PI) / 180, c);
};
A.EqCoord = function (a, b) {
    if (isNaN(a) || isNaN(b))
        throw Error('Invalid EqCoord object: (' + a + ', ' + b + ')');
    this.ra = a;
    this.dec = b;
};
A.EqCoord.prototype = {
    toString: function () {
        return 'ra:' + A.Math.formatNum((180 * this.ra) / Math.PI) + ', dec:' + A.Math.formatNum((180 * this.dec) / Math.PI);
    },
};
A.HzCoord = function (a, b) {
    if (isNaN(a) || isNaN(b))
        throw Error('Invalid HzCoord object: (' + a + ', ' + b + ')');
    this.az = a;
    this.alt = b;
};
A.HzCoord.prototype = {
    toString: function () {
        return 'azi:' + A.Math.formatNum((180 * this.az) / Math.PI) + ', alt:' + A.Math.formatNum((180 * this.alt) / Math.PI);
    },
};
A.Coord = {
    dmsToDeg: function (a, b, c, d) {
        d = (60 * (60 * b + c) + d) / 3600;
        return a ? -d : d;
    },
    calcAngle: function (a, b, c, d) {
        return (A.Coord.dmsToDeg(a, b, c, d) * Math.PI) / 180;
    },
    calcRA: function (a, b, c) {
        return ((A.Coord.dmsToDeg(!1, a, b, c) % 24) * 15 * Math.PI) / 180;
    },
    secondsToHMSStr: function (a) {
        var b = Math.floor(a / 86400);
        a = A.Math.pMod(a, 86400);
        var c = Math.floor(a / 3600) % 24, d = Math.floor(a / 60) % 60;
        a = Math.floor(a % 60);
        return (0 !== b ? b + 'd ' : '') + (10 > c ? '0' : '') + c + ':' + (10 > d ? '0' : '') + d + ':' + (10 > a ? '0' : '') + a;
    },
    secondsToHMStr: function (a) {
        var b = Math.floor(a / 86400);
        a = A.Math.pMod(a, 86400);
        var c = Math.floor(a / 3600) % 24;
        a = Math.floor(a / 60) % 60;
        return (0 !== b ? b + 'd ' : '') + (10 > c ? '0' : '') + c + ':' + (10 > a ? '0' : '') + a;
    },
    eqToEcl: function (a, b) {
        var c = Math.sin(a.ra), d = Math.sin(a.dec), e = Math.cos(a.dec), f = Math.sin(b);
        b = Math.cos(b);
        return new A.EclCoord(Math.atan2(c * b + (d / e) * f, Math.cos(a.ra)), Math.asin(d * b - e * f * c));
    },
    eclToEq: function (a, b) {
        var c = Math.sin(a.lat), d = Math.sin(a.lng), e = Math.cos(a.lng), f = Math.sin(b);
        b = Math.cos(b);
        let a2 = Math.atan2(c * b - (d / e) * f, Math.cos(a.lat));
        0 > a2 && (a2 += 2 * Math.PI);
        return new A.EqCoord(a2, Math.asin(d * b + e * f * c));
    },
    eqToHz: function (a, b, c) {
        c = c - b.lng - a.ra;
        var d = Math.cos(c), e = Math.sin(b.lat);
        b = Math.cos(b.lat);
        var f = Math.sin(a.dec);
        a = Math.cos(a.dec);
        return new A.HzCoord(Math.atan2(Math.sin(c), d * e - (f / a) * b), Math.asin(e * f + b * a * d));
    },
};
A.DeltaT = {
    jdToJde: function (a, b) {
        b || (b = A.DeltaT.estimate(a));
        return a + b / 86400;
    },
    jdeToJd: function (a, b) {
        b || (b = A.DeltaT.estimate(a));
        return a - b / 86400;
    },
    decimalYear: function (a) {
        a = A.JulianDay.jdToCalendar(a);
        return a.y + (a.m - 0.5) / 12;
    },
    estimate: function (a) {
        var b = A.DeltaT.decimalYear(a);
        a = Math.pow;
        return -500 > b
            ? -20 + 32 * a((b - 1820) / 100, 2)
            : 500 > b
                ? ((b /= 100), 10583.6 - 1014.41 * b + 33.78311 * a(b, 2) - 5.952053 * a(b, 3) - 0.1798452 * a(b, 4) + 0.022174192 * a(b, 5) + 0.0090316521 * a(b, 6))
                : 1600 > b
                    ? ((b = (b - 1e3) / 100), 1574.2 - 556.01 * b + 71.23472 * a(b, 2) + 0.319781 * a(b, 3) - 0.8503463 * a(b, 4) - 0.005050998 * a(b, 5) + 0.0083572073 * a(b, 6))
                    : 1700 > b
                        ? ((b -= 1600), 120 - 0.9808 * b - 0.01532 * a(b, 2) + a(b, 3) / 7129)
                        : 1800 > b
                            ? ((b -= 1700), 8.83 + 0.1603 * b - 0.0059285 * a(b, 2) + 1.3336e-4 * a(b, 3) - a(b, 4) / 1174e3)
                            : 1860 > b
                                ? ((b -= 1800), 13.72 - 0.332447 * b + 0.0068612 * a(b, 2) + 0.0041116 * a(b, 3) - 3.7436e-4 * a(b, 4) + 1.21272e-5 * a(b, 5) - 1.699e-7 * a(b, 6) + 8.75e-10 * a(b, 7))
                                : 1900 > b
                                    ? ((b -= 1860), 7.62 + 0.5737 * b - 0.251754 * a(b, 2) + 0.01680668 * a(b, 3) - 4.473624e-4 * a(b, 4) + a(b, 5) / 233174)
                                    : 1920 > b
                                        ? ((b -= 1900), -2.79 + 1.494119 * b - 0.0598939 * a(b, 2) + 0.0061966 * a(b, 3) - 1.97e-4 * a(b, 4))
                                        : 1941 > b
                                            ? ((b -= 1920), 21.2 + 0.84493 * b - 0.0761 * a(b, 2) + 0.0020936 * a(b, 3))
                                            : 1961 > b
                                                ? ((b -= 1950), 29.07 + 0.407 * b - a(b, 2) / 233 + a(b, 3) / 2547)
                                                : 1986 > b
                                                    ? ((b -= 1975), 45.45 + 1.067 * b - a(b, 2) / 260 - a(b, 3) / 718)
                                                    : 2005 > b
                                                        ? ((b -= 2e3), 63.86 + 0.3345 * b - 0.060374 * a(b, 2) + 0.0017275 * a(b, 3) + 6.51814e-4 * a(b, 4) + 2.373599e-5 * a(b, 5))
                                                        : 2050 > b
                                                            ? ((b -= 2e3), 62.92 + 0.32217 * b + 0.005589 * a(b, 2))
                                                            : 2150 > b
                                                                ? -20 + 32 * a((b - 1820) / 100, 2) - 0.5628 * (2150 - b)
                                                                : -20 + 32 * a((b - 1820) / 100, 2);
    },
};
A.Globe = {
    Er: 6378.14,
    Fl: 1 / 298.257,
    parallaxConstants: function (a, b) {
        b || (b = 0);
        var c = 1 - A.Globe.Fl;
        b = (0.001 * b) / A.Globe.Er;
        return { rhoslat: Math.sin(Math.atan(c * Math.tan(a))) * c + b * Math.sin(a), rhoclat: Math.cos(Math.atan(c * Math.tan(a))) + b * Math.cos(a) };
    },
};
A.Interp = {
    newLen3: function (a, b, c) {
        if (3 != c.length)
            throw 'Error not 3';
        if (b === a)
            throw 'Error no x range';
        var d = c[1] - c[0], e = c[2] - c[1];
        return { x1: a, x3: b, y: c, a: d, b: e, c: e - d, abSum: d + e, xSum: b + a, xDiff: b - a };
    },
    interpolateX: function (a, b) {
        return A.Interp.interpolateN(a, (2 * b - a.xSum) / a.xDiff);
    },
    interpolateN: function (a, b) {
        return a.y[1] + 0.5 * b * (a.abSum + b * a.c);
    },
};
A.JulianDay = function (a, b) {
    a instanceof Date && (a = A.JulianDay.dateToJD(a));
    this.jd = a;
    this.deltaT = b ? b : A.DeltaT.estimate(this.jd);
    this.jde = A.DeltaT.jdToJde(this.jd, this.deltaT);
};
A.JulianDay.prototype = {
    toCalendar: function () {
        return A.JulianDay.jdToCalendar(this.jd);
    },
    toDate: function () {
        return A.JulianDay.jdToDate(this.jd);
    },
    jdJ2000Century: function () {
        return (this.jd - A.J2000) / A.JulianCentury;
    },
    jdeJ2000Century: function () {
        return (this.jde - A.J2000) / A.JulianCentury;
    },
    startOfDay: function () {
        return new A.JulianDay(Math.floor(this.jde - 0.5) + 0.5, this.deltaT);
    },
};
A.JulianDay.gregorianTimeStart = Date.UTC(1582, 9, 4);
A.JulianDay.jdFromGregorian = function (a, b, c) {
    return new A.JulianDay(A.JulianDay.jdFromGregorian(a, b, c));
};
A.JulianDay.jdFromJulian = function (a, b, c) {
    return new A.JulianDay(A.JulianDay.calendarJulianToJD(a, b, c));
};
A.JulianDay.jdFromJDE = function (a) {
    var b = A.DeltaT.estimate(a);
    a = A.DeltaT.jdeToJd(a, b);
    return new A.JulianDay(a, b);
};
A.JulianDay.dateToJD = function (a) {
    var b = a.getUTCDate() + A.JulianDay.secondsFromHMS(a.getUTCHours(), a.getUTCMinutes(), a.getUTCSeconds()) / 86400;
    return a.getTime() < A.JulianDay.gregorianTimeStart
        ? A.JulianDay.calendarJulianToJD(a.getUTCFullYear(), a.getUTCMonth() + 1, b)
        : A.JulianDay.calendarGregorianToJD(a.getUTCFullYear(), a.getUTCMonth() + 1, b);
};
A.JulianDay.calendarGregorianToJD = function (a, b, c) {
    if (1 === b || 2 === b)
        a--, (b += 12);
    var d = Math.floor(a / 100);
    return Math.floor((36525 * (a + 4716)) / 100) + Math.floor((306 * (b + 1)) / 10) + (2 - d + Math.floor(d / 4)) + c - 1524.5;
};
A.JulianDay.calendarJulianToJD = function (a, b, c) {
    if (1 === b || 2 === b)
        a--, (b += 12);
    return Math.floor((36525 * (a + 4716)) / 100) + Math.floor((306 * (b + 1)) / 10) + c - 1524.5;
};
A.JulianDay.secondsFromHMS = function (a, b, c) {
    return 3600 * a + 60 * b + c;
};
A.JulianDay.jdToDate = function (a) {
    var b = A.JulianDay.jdToCalendar(a);
    a = A.Math.modF(a + 0.5)[1];
    a = Math.round(86400 * a);
    return new Date(Date.UTC(b.y, b.m - 1, Math.floor(b.d), Math.floor(a / 3600) % 24, Math.floor(a / 60) % 60, Math.floor(a % 60)));
};
A.JulianDay.jdToCalendar = function (a) {
    a = A.Math.modF(a + 0.5);
    var b = a[0], c = b;
    2299151 <= b && ((c = Math.floor((100 * b - 186721625) / 3652425)), (c = b + 1 + c - Math.floor(c / 4)));
    var d = c + 1524;
    b = Math.floor((100 * d - 12210) / 36525);
    var e = Math.floor((36525 * b) / 100);
    c = Math.floor((1e4 * (d - e)) / 306001);
    a = d - e - Math.floor((306001 * c) / 1e4) + a[1];
    c = 14 === c || 15 === c ? c - 13 : c - 1;
    return { y: 1 === c || 2 === c ? Math.floor(b) - 4715 : Math.floor(b) - 4716, m: c, d: a };
};
A.JulianDay.leapYearGregorian = function (a) {
    return (0 === a % 4 && 0 !== a % 100) || 0 === a % 400;
};
A.JulianDay.dayOfYear = function (a, b, c, d) {
    a = 2;
    d && a--;
    return A.JulianDay._wholeMonths(b, a) + c;
};
A.JulianDay._wholeMonths = function (a, b) {
    return Math.round((275 * a) / 9 - ((a + 9) / 12) * b - 30);
};
A.Math = {
    pMod: function (a, b) {
        a %= b;
        0 > a && (a += b);
        return a;
    },
    modF: function (a) {
        return 0 > a ? ((a = -a), [-Math.floor(a), -(a % 1)]) : [Math.floor(a), a % 1];
    },
    horner: function (a, b) {
        var c = b.length - 1;
        if (0 >= c)
            throw 'empty array not supported';
        for (var d = b[c]; 0 < c;)
            c--, (d = d * a + b[c]);
        return d;
    },
    formatNum: function (a, b) {
        b = Math.pow(10, b | 4);
        return Math.round(a * b) / b;
    },
};
A.Moon = {
    parallax: function (a) {
        return Math.asin(6378.14 / a);
    },
    apparentEquatorial: function (a) {
        var b = A.Moon.geocentricPosition(a), c = A.Nutation.nutation(a);
        a = A.Nutation.meanObliquityLaskar(a) + c.deltaobliquity;
        return { eq: A.Coord.eclToEq(new A.EclCoord(b.lng + c.deltalng, b.lat), a), delta: b.delta };
    },
    apparentTopocentric: function (a, b, c) {
        var d = A.Moon.apparentEquatorial(a), e = A.Globe.parallaxConstants(b.lat, b.h), f = A.Moon.parallax(d.delta);
        c || (c = A.Sidereal.apparentInRa(a));
        return { eq: A.Parallax.topocentric(d.eq, f, e.rhoslat, e.rhoclat, b.lng, c), delta: d.delta };
    },
    topocentricPosition: function (a, b, c) {
        var d = A.Sidereal.apparentInRa(a);
        a = A.Moon.apparentTopocentric(a, b, d);
        var e = A.Coord.eqToHz(a.eq, b, d);
        !0 === c && (e.alt += A.Refraction.bennett2(e.alt));
        b = A.Moon.parallacticAngle(b.lat, d - (b.lng + a.eq.ra), a.eq.dec);
        return { hz: e, eq: a.eq, delta: a.delta, q: b };
    },
    approxTransit: function (a, b) {
        a = a.startOfDay();
        return A.Rise.approxTransit(b, A.Sidereal.apparent0UT(a), A.Moon.apparentTopocentric(a, b).eq);
    },
    approxTimes: function (a, b) {
        a = a.startOfDay();
        var c = A.Moon.apparentTopocentric(a, b), d = A.Moon.parallax(c.delta);
        d = A.Rise.stdh0Lunar(d);
        a = A.Sidereal.apparent0UT(a);
        return A.Rise.approxTimes(b, d, a, c.eq);
    },
    times: function (a, b) {
        a = a.startOfDay();
        var c = A.Moon.apparentTopocentric(new A.JulianDay(a.jd - 1, a.deltaT), b), d = A.Moon.apparentTopocentric(a, b), e = A.Moon.apparentTopocentric(new A.JulianDay(a.jd + 1, a.deltaT), b), f = A.Moon.parallax(d.delta);
        f = A.Rise.stdh0Lunar(f);
        var g = A.Sidereal.apparent0UT(a);
        return A.Rise.times(b, a.deltaT, f, g, [c.eq, d.eq, e.eq]);
    },
    parallacticAngle: function (a, b, c) {
        return Math.atan2(Math.sin(b), Math.tan(a) * Math.cos(c) - Math.sin(c) * Math.cos(b));
    },
    geocentricPosition: function (a) {
        var b = Math.PI / 180, c = a.jdeJ2000Century();
        a = A.Math.pMod(A.Math.horner(c, [218.3164477 * b, 481267.88123421 * b, -0.0015786 * b, b / 538841, -b / 65194e3]), 2 * Math.PI);
        var d = A.Math.pMod(A.Math.horner(c, [297.8501921 * b, 445267.1114034 * b, -0.0018819 * b, b / 545868, -b / 113065e3]), 2 * Math.PI), e = A.Math.pMod(A.Math.horner(c, [357.5291092 * b, 35999.0502909 * b, -1.535e-4 * b, b / 2449e4]), 2 * Math.PI), f = A.Math.pMod(A.Math.horner(c, [134.9633964 * b, 477198.8675055 * b, 0.0087414 * b, b / 69699, -b / 14712e3]), 2 * Math.PI), g = A.Math.pMod(A.Math.horner(c, [93.272095 * b, 483202.0175233 * b, -0.0036539 * b, -b / 3526e3, b / 86331e4]), 2 * Math.PI), l = 119.75 * b + 131.849 * b * c, m = 53.09 * b + 479264.29 * b * c, h = 313.45 * b + 481266.484 * b * c;
        c = A.Math.horner(c, [1, -0.002516, -7.4e-6]);
        var p = c * c;
        m = 3958 * Math.sin(l) + 1962 * Math.sin(a - g) + 318 * Math.sin(m);
        var n = 0;
        l = -2235 * Math.sin(a) + 382 * Math.sin(h) + 175 * Math.sin(l - g) + 175 * Math.sin(l + g) + 127 * Math.sin(a - f) - 115 * Math.sin(a + f);
        for (h = 0; h < A.Moon.ta.length; h++) {
            var k = A.Moon.ta[h];
            var r = d * k[0] + e * k[1] + f * k[2] + g * k[3], q = Math.sin(r);
            r = Math.cos(r);
            switch (k[1]) {
                case 0:
                    m += k[4] * q;
                    n += k[5] * r;
                    break;
                case 1:
                case -1:
                    m += k[4] * q * c;
                    n += k[5] * r * c;
                    break;
                case 2:
                case -2:
                    m += k[4] * q * p;
                    n += k[5] * r * p;
                    break;
                default:
                    throw 'error';
            }
        }
        for (h = 0; h < A.Moon.tb.length; h++)
            switch (((k = A.Moon.tb[h]), (q = Math.sin(d * k[0] + e * k[1] + f * k[2] + g * k[3])), k[1])) {
                case 0:
                    l += k[4] * q;
                    break;
                case 1:
                case -1:
                    l += k[4] * q * c;
                    break;
                case 2:
                case -2:
                    l += k[4] * q * p;
                    break;
                default:
                    throw 'error';
            }
        return { lng: A.Math.pMod(a, 2 * Math.PI) + 1e-6 * m * b, lat: 1e-6 * l * b, delta: 385000.56 + 0.001 * n };
    },
    ta: [
        [0, 0, 1, 0, 6288774, -20905355],
        [2, 0, -1, 0, 1274027, -3699111],
        [2, 0, 0, 0, 658314, -2955968],
        [0, 0, 2, 0, 213618, -569925],
        [0, 1, 0, 0, -185116, 48888],
        [0, 0, 0, 2, -114332, -3149],
        [2, 0, -2, 0, 58793, 246158],
        [2, -1, -1, 0, 57066, -152138],
        [2, 0, 1, 0, 53322, -170733],
        [2, -1, 0, 0, 45758, -204586],
        [0, 1, -1, 0, -40923, -129620],
        [1, 0, 0, 0, -34720, 108743],
        [0, 1, 1, 0, -30383, 104755],
        [2, 0, 0, -2, 15327, 10321],
        [0, 0, 1, 2, -12528, 0],
        [0, 0, 1, -2, 10980, 79661],
        [4, 0, -1, 0, 10675, -34782],
        [0, 0, 3, 0, 10034, -23210],
        [4, 0, -2, 0, 8548, -21636],
        [2, 1, -1, 0, -7888, 24208],
        [2, 1, 0, 0, -6766, 30824],
        [1, 0, -1, 0, -5163, -8379],
        [1, 1, 0, 0, 4987, -16675],
        [2, -1, 1, 0, 4036, -12831],
        [2, 0, 2, 0, 3994, -10445],
        [4, 0, 0, 0, 3861, -11650],
        [2, 0, -3, 0, 3665, 14403],
        [0, 1, -2, 0, -2689, -7003],
        [2, 0, -1, 2, -2602, 0],
        [2, -1, -2, 0, 2390, 10056],
        [1, 0, 1, 0, -2348, 6322],
        [2, -2, 0, 0, 2236, -9884],
        [0, 1, 2, 0, -2120, 5751],
        [0, 2, 0, 0, -2069, 0],
        [2, -2, -1, 0, 2048, -4950],
        [2, 0, 1, -2, -1773, 4130],
        [2, 0, 0, 2, -1595, 0],
        [4, -1, -1, 0, 1215, -3958],
        [0, 0, 2, 2, -1110, 0],
        [3, 0, -1, 0, -892, 3258],
        [2, 1, 1, 0, -810, 2616],
        [4, -1, -2, 0, 759, -1897],
        [0, 2, -1, 0, -713, -2117],
        [2, 2, -1, 0, -700, 2354],
        [2, 1, -2, 0, 691, 0],
        [2, -1, 0, -2, 596, 0],
        [4, 0, 1, 0, 549, -1423],
        [0, 0, 4, 0, 537, -1117],
        [4, -1, 0, 0, 520, -1571],
        [1, 0, -2, 0, -487, -1739],
        [2, 1, 0, -2, -399, 0],
        [0, 0, 2, -2, -381, -4421],
        [1, 1, 1, 0, 351, 0],
        [3, 0, -2, 0, -340, 0],
        [4, 0, -3, 0, 330, 0],
        [2, -1, 2, 0, 327, 0],
        [0, 2, 1, 0, -323, 1165],
        [1, 1, -1, 0, 299, 0],
        [2, 0, 3, 0, 294, 0],
        [2, 0, -1, -2, 0, 8752],
    ],
    tb: [
        [0, 0, 0, 1, 5128122],
        [0, 0, 1, 1, 280602],
        [0, 0, 1, -1, 277693],
        [2, 0, 0, -1, 173237],
        [2, 0, -1, 1, 55413],
        [2, 0, -1, -1, 46271],
        [2, 0, 0, 1, 32573],
        [0, 0, 2, 1, 17198],
        [2, 0, 1, -1, 9266],
        [0, 0, 2, -1, 8822],
        [2, -1, 0, -1, 8216],
        [2, 0, -2, -1, 4324],
        [2, 0, 1, 1, 4200],
        [2, 1, 0, -1, -3359],
        [2, -1, -1, 1, 2463],
        [2, -1, 0, 1, 2211],
        [2, -1, -1, -1, 2065],
        [0, 1, -1, -1, -1870],
        [4, 0, -1, -1, 1828],
        [0, 1, 0, 1, -1794],
        [0, 0, 0, 3, -1749],
        [0, 1, -1, 1, -1565],
        [1, 0, 0, 1, -1491],
        [0, 1, 1, 1, -1475],
        [0, 1, 1, -1, -1410],
        [0, 1, 0, -1, -1344],
        [1, 0, 0, -1, -1335],
        [0, 0, 3, 1, 1107],
        [4, 0, 0, -1, 1021],
        [4, 0, -1, 1, 833],
        [0, 0, 1, -3, 777],
        [4, 0, -2, 1, 671],
        [2, 0, 0, -3, 607],
        [2, 0, 2, -1, 596],
        [2, -1, 1, -1, 491],
        [2, 0, -2, 1, -451],
        [0, 0, 3, -1, 439],
        [2, 0, 2, 1, 422],
        [2, 0, -3, -1, 421],
        [2, 1, -1, 1, -366],
        [2, 1, 0, 1, -351],
        [4, 0, 0, 1, 331],
        [2, -1, 1, 1, 315],
        [2, -2, 0, -1, 302],
        [0, 0, 1, 3, -283],
        [2, 1, 1, -1, -229],
        [1, 1, 0, -1, 223],
        [1, 1, 0, 1, 223],
        [0, 1, -2, -1, -220],
        [2, 1, -1, -1, -220],
        [1, 0, 1, 1, -185],
        [2, -1, -2, -1, 181],
        [0, 1, 2, 1, -177],
        [4, 0, -2, -1, 176],
        [4, -1, -1, -1, 166],
        [1, 0, 1, -1, -164],
        [4, 0, 1, -1, 132],
        [1, 0, -1, -1, -119],
        [4, -1, 0, -1, 115],
        [2, -2, 0, 1, 107],
    ],
};
A.MoonIllum = {
    phaseAngleEq: function (a, b, c, d) {
        a = A.MoonIllum._coselong(a, c);
        return Math.atan2(d * Math.sin(Math.acos(a)), b - d * a);
    },
    phaseAngleEq2: function (a, b) {
        return Math.acos(-A.MoonIllum._coselong(a, b));
    },
    illuminated: function (a) {
        return (1 + Math.cos(a)) / 2;
    },
    positionAngle: function (a, b) {
        var c = Math.cos(b.dec);
        return Math.atan2(c * Math.sin(b.ra - a.ra), Math.sin(b.dec) * Math.cos(a.dec) - c * Math.sin(a.dec) * Math.cos(b.ra - a.ra));
    },
    _coselong: function (a, b) {
        return Math.sin(b.dec) * Math.sin(a.dec) + Math.cos(b.dec) * Math.cos(a.dec) * Math.cos(b.ra - a.ra);
    },
};
A.Nutation = {
    nutation: function (a) {
        a = a.jdeJ2000Century();
        for (var b = (A.Math.horner(a, [297.85036, 445267.11148, -0.0019142, 1 / 189474]) * Math.PI) / 180, c = (A.Math.horner(a, [357.52772, 35999.05034, -1.603e-4, -1 / 3e5]) * Math.PI) / 180, d = (A.Math.horner(a, [134.96298, 477198.867398, 0.0086972, 1 / 5620]) * Math.PI) / 180, e = (A.Math.horner(a, [93.27191, 483202.017538, -0.0036825, 1 / 327270]) * Math.PI) / 180, f = (A.Math.horner(a, [125.04452, -1934.136261, 0.0020708, 1 / 45e4]) * Math.PI) / 180, g = 0, l = 0, m = A.Nutation.table22A.length - 1; 0 <= m; m--) {
            var h = A.Nutation.table22A[m], p = h[0] * b + h[1] * c + h[2] * d + h[3] * e + h[4] * f, n = Math.cos(p);
            g += Math.sin(p) * (h[5] + h[6] * a);
            l += n * (h[7] + h[8] * a);
        }
        return { deltalng: (1e-4 / 3600) * g * (Math.PI / 180), deltaobliquity: (1e-4 / 3600) * l * (Math.PI / 180) };
    },
    nutationInRA: function (a) {
        var b = A.Nutation.meanObliquityLaskar(a);
        a = A.Nutation.nutation(a);
        return a.deltalng * Math.cos(b + a.deltaobliquity);
    },
    trueObliquity: function (a) {
        var b = A.Nutation.meanObliquityLaskar(a);
        a = A.Nutation.nutation(a);
        return b + a.deltaobliquity;
    },
    meanObliquity: function (a) {
        return A.Math.horner(a.jdeJ2000Century(), [
            (84381.448 / 3600) * (Math.PI / 180),
            (-46.815 / 3600) * (Math.PI / 180),
            (-5.9e-4 / 3600) * (Math.PI / 180),
            (0.001813 / 3600) * (Math.PI / 180),
        ]);
    },
    meanObliquityLaskar: function (a) {
        return A.Math.horner(0.01 * a.jdeJ2000Century(), [
            (84381.448 / 3600) * (Math.PI / 180),
            (-4680.93 / 3600) * (Math.PI / 180),
            (-1.55 / 3600) * (Math.PI / 180),
            (1999.25 / 3600) * (Math.PI / 180),
            (-51.38 / 3600) * (Math.PI / 180),
            (-249.67 / 3600) * (Math.PI / 180),
            (-39.05 / 3600) * (Math.PI / 180),
            (7.12 / 3600) * (Math.PI / 180),
            (27.87 / 3600) * (Math.PI / 180),
            (5.79 / 3600) * (Math.PI / 180),
            (2.45 / 3600) * (Math.PI / 180),
        ]);
    },
    table22A: [
        [0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9],
        [-2, 0, 0, 2, 2, -13187, -1.6, 5736, -3.1],
        [0, 0, 0, 2, 2, -2274, -0.2, 977, -0.5],
        [0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5],
        [0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1],
        [0, 0, 1, 0, 0, 712, 0.1, -7, 0],
        [-2, 1, 0, 2, 2, -517, 1.2, 224, -0.6],
        [0, 0, 0, 2, 1, -386, -0.4, 200, 0],
        [0, 0, 1, 2, 2, -301, 0, 129, -0.1],
        [-2, -1, 0, 2, 2, 217, -0.5, -95, 0.3],
        [-2, 0, 1, 0, 0, -158, 0, 0, 0],
        [-2, 0, 0, 2, 1, 129, 0.1, -70, 0],
        [0, 0, -1, 2, 2, 123, 0, -53, 0],
        [2, 0, 0, 0, 0, 63, 0, 0, 0],
        [0, 0, 1, 0, 1, 63, 0.1, -33, 0],
        [2, 0, -1, 2, 2, -59, 0, 26, 0],
        [0, 0, -1, 0, 1, -58, -0.1, 32, 0],
        [0, 0, 1, 2, 1, -51, 0, 27, 0],
        [-2, 0, 2, 0, 0, 48, 0, 0, 0],
        [0, 0, -2, 2, 1, 46, 0, -24, 0],
        [2, 0, 0, 2, 2, -38, 0, 16, 0],
        [0, 0, 2, 2, 2, -31, 0, 13, 0],
        [0, 0, 2, 0, 0, 29, 0, 0, 0],
        [-2, 0, 1, 2, 2, 29, 0, -12, 0],
        [0, 0, 0, 2, 0, 26, 0, 0, 0],
        [-2, 0, 0, 2, 0, -22, 0, 0, 0],
        [0, 0, -1, 2, 1, 21, 0, -10, 0],
        [0, 2, 0, 0, 0, 17, -0.1, 0, 0],
        [2, 0, -1, 0, 1, 16, 0, -8, 0],
        [-2, 2, 0, 2, 2, -16, 0.1, 7, 0],
        [0, 1, 0, 0, 1, -15, 0, 9, 0],
        [-2, 0, 1, 0, 1, -13, 0, 7, 0],
        [0, -1, 0, 0, 1, -12, 0, 6, 0],
        [0, 0, 2, -2, 0, 11, 0, 0, 0],
        [2, 0, -1, 2, 1, -10, 0, 5, 0],
        [2, 0, 1, 2, 2, -8, 0, 3, 0],
        [0, 1, 0, 2, 2, 7, 0, -3, 0],
        [-2, 1, 1, 0, 0, -7, 0, 0, 0],
        [0, -1, 0, 2, 2, -7, 0, 3, 0],
        [2, 0, 0, 2, 1, -7, 0, 3, 0],
        [2, 0, 1, 0, 0, 6, 0, 0, 0],
        [-2, 0, 2, 2, 2, 6, 0, -3, 0],
        [-2, 0, 1, 2, 1, 6, 0, -3, 0],
        [2, 0, -2, 0, 1, -6, 0, 3, 0],
        [2, 0, 0, 0, 1, -6, 0, 3, 0],
        [0, -1, 1, 0, 0, 5, 0, 0, 0],
        [-2, -1, 0, 2, 1, -5, 0, 3, 0],
        [-2, 0, 0, 0, 1, -5, 0, 3, 0],
        [0, 0, 2, 2, 1, -5, 0, 3, 0],
        [-2, 0, 2, 0, 1, 4, 0, 0, 0],
        [-2, 1, 0, 2, 1, 4, 0, 0, 0],
        [0, 0, 1, -2, 0, 4, 0, 0, 0],
        [-1, 0, 1, 0, 0, -4, 0, 0, 0],
        [-2, 1, 0, 0, 0, -4, 0, 0, 0],
        [1, 0, 0, 0, 0, -4, 0, 0, 0],
        [0, 0, 1, 2, 0, 3, 0, 0, 0],
        [0, 0, -2, 2, 2, -3, 0, 0, 0],
        [-1, -1, 1, 0, 0, -3, 0, 0, 0],
        [0, 1, 1, 0, 0, -3, 0, 0, 0],
        [0, -1, 1, 2, 2, -3, 0, 0, 0],
        [2, -1, -1, 2, 2, -3, 0, 0, 0],
        [0, 0, 3, 2, 2, -3, 0, 0, 0],
        [2, -1, 0, 2, 2, -3, 0, 0, 0],
    ],
};
A.Parallax = {
    earthsunParallax: ((8.794 / 60 / 60) * Math.PI) / 180,
    horizontal: function (a) {
        return ((8.794 / 60 / 60) * Math.PI) / 180 / a;
    },
    topocentric: function (a, b, c, d, e, f) {
        e = A.Math.pMod(f - e - a.ra, 2 * Math.PI);
        b = Math.sin(b);
        f = Math.cos(e);
        var g = Math.cos(a.dec);
        e = Math.atan2(-d * b * Math.sin(e), g - d * b * f);
        return new A.EqCoord(a.ra + e, Math.atan2((Math.sin(a.dec) - c * b) * Math.cos(e), g - d * b * f));
    },
    topocentric2: function (a, b, c, d, e, f) {
        e = A.Math.pMod(f - e - a.ra, 2 * Math.PI);
        f = Math.cos(a.dec);
        return new A.EqCoord(a.ra + (-b * d * Math.sin(e)) / f, a.dec + -b * (c * f - d * Math.cos(e) * Math.sin(a.dec)));
    },
};
A.Refraction = {
    bennett: function (a) {
        0 > a && (a = 0);
        var b = Math.PI / 180;
        return b / 60 / Math.tan(a + (7.31 * b * b) / (a + 4.4 * b));
    },
    bennett2: function (a) {
        var b = Math.PI / 180, c = 60 / b, d = 0.06 / c;
        c = 14.7 * c * b;
        b *= 13;
        a = A.Refraction.bennett(a);
        return a - d * Math.sin(c * a + b);
    },
    saemundsson: function (a) {
        var b = Math.PI / 180;
        return (1.02 * b) / 60 / Math.tan(a + (10.3 * b * b) / (a + 5.11 * b));
    },
};
A.Rise = {
    meanRefraction: (0.5667 * Math.PI) / 180,
    stdh0Stellar: (-0.5667 * Math.PI) / 180,
    stdh0Solar: (-0.8333 * Math.PI) / 180,
    stdh0LunarMean: (0.125 * Math.PI) / 180,
    stdh0Lunar: function (a) {
        return 0.7275 * a - A.Rise.meanRefraction;
    },
    circumpolar: function (a, b, c) {
        a = (Math.sin(b) - Math.sin(a) * Math.sin(c)) / (Math.cos(a) * Math.cos(c));
        return -1 > a || 1 < a ? null : a;
    },
    approxTransit: function (a, b, c) {
        return (43200 * (c.ra + a.lng)) / Math.PI - b;
    },
    approxTimes: function (a, b, c, d) {
        b = A.Rise.circumpolar(a.lat, b, d.dec);
        if (!b)
            return null;
        b = (43200 * Math.acos(b)) / Math.PI;
        a = (43200 * (d.ra + a.lng)) / Math.PI - c;
        return {
            transit: A.Math.pMod(a, 86400),
            transitd: Math.floor(a / 86400),
            rise: A.Math.pMod(a - b, 86400),
            rised: Math.floor((a - b) / 86400),
            set: A.Math.pMod(a + b, 86400),
            setd: Math.floor((a + b) / 86400),
        };
    },
    times: function (a, b, c, d, e) {
        function f(e) {
            var f = A.Math.pMod(d + (360.985647 * e) / 360, 86400), g = e + b, h = A.Interp.interpolateX(l, g);
            g = A.Interp.interpolateX(m, g);
            f = (f * Math.PI) / 43200 - (a.lng + h);
            h = Math.cos(g);
            return A.Math.pMod(e + (((p * Math.sin(g) + n * h * Math.cos(f) - c) / (h * n * Math.sin(f))) * 43200) / Math.PI, 86400);
        }
        var g = A.Rise.approxTimes(a, c, d, e[1]);
        if (!g)
            return null;
        var l = A.Interp.newLen3(-86400, 86400, [e[0].ra, e[1].ra, e[2].ra]), m = A.Interp.newLen3(-86400, 86400, [e[0].dec, e[1].dec, e[2].dec]);
        e = d + (360.985647 * g.transit) / 360;
        var h = A.Interp.interpolateX(l, g.transit + b);
        g.transit = A.Math.pMod(g.transit - (e - (43200 * (a.lng + h)) / Math.PI), 86400);
        var p = Math.sin(a.lat), n = Math.cos(a.lat);
        g.rise = f(g.rise);
        g.set = f(g.set);
        return g;
    },
};
A.Sidereal = {
    iau82: [24110.54841, 8640184.812866, 0.093104, 6.2e-6],
    jdToCFrac: function (a) {
        a = A.Math.modF(a.jd + 0.5);
        return [new A.JulianDay(a[0] - 0.5).jdJ2000Century(), a[1]];
    },
    mean: function (a) {
        return A.Math.pMod(A.Sidereal._mean(a), 86400);
    },
    _mean: function (a) {
        a = A.Sidereal._mean0UT(a);
        return a.s + 86636.55536784 * a.f;
    },
    _meanInRA: function (a) {
        a = A.Sidereal._mean0UT(a);
        return (a.s * Math.PI) / 43200 + 2.0054758187 * a.f * Math.PI;
    },
    mean0UT: function (a) {
        a = A.Sidereal._mean0UT(a);
        return A.Math.pMod(a.s, 86400);
    },
    _mean0UT: function (a) {
        a = A.Sidereal.jdToCFrac(a);
        return { s: A.Math.horner(a[0], A.Sidereal.iau82), f: a[1] };
    },
    apparentInRa: function (a) {
        var b = A.Sidereal._meanInRA(a);
        a = A.Nutation.nutationInRA(a);
        return A.Math.pMod(b + a, 2 * Math.PI);
    },
    apparent: function (a) {
        var b = A.Sidereal._mean(a);
        a = (648e3 * A.Nutation.nutationInRA(a)) / Math.PI / 15;
        return A.Math.pMod(b + a, 86400);
    },
    apparentLocal: function (a, b) {
        a = A.Sidereal.apparent(a);
        return A.Math.pMod(a - (43200 * b) / Math.PI, 86400);
    },
    apparent0UT: function (a) {
        var b = A.Math.modF(a.jd + 0.5);
        a = A.Math.modF(a.jde + 0.5);
        b = A.Math.horner((b[0] - 0.5 - A.J2000) / 36525, A.Sidereal.iau82) + 86636.55536784 * b[1];
        a = (648e3 * A.Nutation.nutationInRA(new A.JulianDay(a[0]))) / Math.PI / 15;
        return A.Math.pMod(b + a, 86400);
    },
};
A.Solar = {
    earthsunDelta: 149597870,
    apparentEquatorial: function (a) {
        var b = a.jdJ2000Century(), c = A.Solar.node(b);
        b = A.Solar.apparentLongitude(b, c);
        a = A.Nutation.meanObliquityLaskar(a) + ((0.00256 * Math.PI) / 180) * Math.cos(c);
        c = Math.sin(b);
        return new A.EqCoord(Math.atan2(Math.cos(a) * c, Math.cos(b)), Math.asin(Math.sin(a) * c));
    },
    apparentTopocentric: function (a, b, c) {
        var d = A.Solar.apparentEquatorial(a), e = A.Globe.parallaxConstants(b.lat, b.h);
        c || (c = A.Sidereal.apparentInRa(a));
        return A.Parallax.topocentric2(d, A.Parallax.earthsunParallax, e.rhoslat, e.rhoclat, b.lng, c);
    },
    topocentricPosition: function (a, b, c) {
        var d = A.Sidereal.apparentInRa(a);
        a = A.Solar.apparentTopocentric(a, b, d);
        b = A.Coord.eqToHz(a, b, d);
        !0 === c && (b.alt += A.Refraction.bennett2(b.alt));
        return { hz: b, eq: a };
    },
    approxTransit: function (a, b) {
        a = a.startOfDay();
        return A.Rise.approxTransit(b, A.Sidereal.apparent0UT(a), A.Solar.apparentTopocentric(a, b));
    },
    approxTimes: function (a, b) {
        var c = a.startOfDay();
        a = A.Solar.apparentTopocentric(c, b);
        var d = A.Rise.stdh0Solar;
        c = A.Sidereal.apparent0UT(c);
        return A.Rise.approxTimes(b, d, c, a);
    },
    times: function (a, b) {
        a = a.startOfDay();
        var c = A.Solar.apparentTopocentric(new A.JulianDay(a.jd - 1, a.deltaT), b), d = A.Solar.apparentTopocentric(a, b), e = A.Solar.apparentTopocentric(new A.JulianDay(a.jd + 1, a.deltaT), b), f = A.Rise.stdh0Solar, g = A.Sidereal.apparent0UT(a);
        return A.Rise.times(b, a.deltaT, f, g, [c, d, e]);
    },
    meanAnomaly: function (a) {
        return (A.Math.horner(a, [357.52911, 35999.05029, -1.537e-4]) * Math.PI) / 180;
    },
    trueLongitude: function (a) {
        var b = (A.Math.horner(a, [280.46646, 36000.76983, 3.032e-4]) * Math.PI) / 180, c = A.Solar.meanAnomaly(a);
        a = ((A.Math.horner(a, [1.914602, -0.004817, -1.4e-5]) * Math.sin(c) + (0.019993 - 1.01e-4 * a) * Math.sin(2 * c) + 2.89e-4 * Math.sin(3 * c)) * Math.PI) / 180;
        return { s: A.Math.pMod(b + a, 2 * Math.PI), v: A.Math.pMod(c + a, 2 * Math.PI) };
    },
    apparentLongitude: function (a, b) {
        b || (b = A.Solar.node(a));
        return A.Solar.trueLongitude(a).s - (0.00569 * Math.PI) / 180 - ((0.00478 * Math.PI) / 180) * Math.sin(b);
    },
    node: function (a) {
        return ((125.04 - 1934.136 * a) * Math.PI) / 180;
    },
};
A.Solistice = {
    march: function (a) {
        return 1e3 > a ? A.Solistice._eq(a, A.Solistice.mc0) : A.Solistice._eq(a - 2e3, A.Solistice.mc2);
    },
    june: function (a) {
        return 1e3 > a ? A.Solistice._eq(a, A.Solistice.jc0) : A.Solistice._eq(a - 2e3, A.Solistice.jc2);
    },
    september: function (a) {
        return 1e3 > a ? A.Solistice._eq(a, A.Solistice.sc0) : A.Solistice._eq(a - 2e3, A.Solistice.sc2);
    },
    december: function (a) {
        return 1e3 > a ? A.Solistice._eq(a, A.Solistice.dc0) : A.Solistice._eq(a - 2e3, A.Solistice.dc2);
    },
    _eq: function (a, b) {
        a = A.Math.horner(0.001 * a, b);
        b = (a - A.J2000) / A.JulianCentury;
        var c = ((35999.373 * Math.PI) / 180) * b - (2.47 * Math.PI) / 180;
        c = 1 + 0.0334 * Math.cos(c) + 7e-4 * Math.cos(2 * c);
        for (var d = 0, e = this.terms.length - 1; 0 <= e; e--) {
            var f = this.terms[e];
            d += f[0] * Math.cos(((f[1] + f[2] * b) * Math.PI) / 180);
        }
        return a + (1e-5 * d) / c;
    },
    mc0: [1721139.29189, 365242.1374, 0.06134, 0.00111, -7.1e-4],
    jc0: [1721233.25401, 365241.72562, -0.05232, 0.00907, 2.5e-4],
    sc0: [1721325.70455, 365242.49558, -0.11677, -0.00297, 7.4e-4],
    dc0: [1721414.39987, 365242.88257, -0.00769, -0.00933, -6e-5],
    mc2: [2451623.80984, 365242.37404, 0.05169, -0.00411, -5.7e-4],
    jc2: [2451716.56767, 365241.62603, 0.00325, 0.00888, -3e-4],
    sc2: [2451810.21715, 365242.01767, -0.11575, 0.00337, 7.8e-4],
    dc2: [2451900.05952, 365242.74049, -0.06223, -0.00823, 3.2e-4],
    terms: [
        [485, 324.96, 1934.136],
        [203, 337.23, 32964.467],
        [199, 342.08, 20.186],
        [182, 27.85, 445267.112],
        [156, 73.14, 45036.886],
        [136, 171.52, 22518.443],
        [77, 222.54, 65928.934],
        [74, 296.72, 3034.906],
        [70, 243.58, 9037.513],
        [58, 119.81, 33718.147],
        [52, 297.17, 150.678],
        [50, 21.02, 2281.226],
        [45, 247.54, 29929.562],
        [44, 325.15, 31555.956],
        [29, 60.93, 4443.417],
        [18, 155.12, 67555.328],
        [17, 288.79, 4562.452],
        [16, 198.04, 62894.029],
        [14, 199.76, 31436.921],
        [12, 95.39, 14577.848],
        [12, 287.11, 31931.756],
        [12, 320.81, 34777.259],
        [9, 227.73, 1222.114],
        [8, 15.45, 16859.074],
    ],
};



}),
"./src/engine/utils/transforms.ts": 
/*!****************************************!*\
  !*** ./src/engine/utils/transforms.ts ***!
  \****************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  alt2zoom: () => (alt2zoom),
  dateFromJday: () => (dateFromJday),
  dateToLocalInIso: () => (dateToLocalInIso),
  jday: () => (jday),
  lat2pitch: () => (lat2pitch),
  localToZulu: () => (localToZulu),
  lon2yaw: () => (lon2yaw),
  normalizeAngle: () => (normalizeAngle)
});
/* ESM import */var _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @ootk/src/main */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _constants__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ./constants */ "./src/engine/utils/constants.ts");


/**
 * This function normalizes the angle to be between -TAU/2 and TAU/2.
 */
const normalizeAngle = (angle) => {
    let normalizedAngle = angle % _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU;
    if (normalizedAngle > _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU / 2) {
        normalizedAngle -= _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU;
    }
    if (normalizedAngle < -_ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU / 2) {
        normalizedAngle += _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU;
    }
    return normalizedAngle;
};
/**
 * This function converts longitude in degrees to yaw in radians for the camera.
 */
const lon2yaw = (lon, selectedDate) => {
    const realTime = new Date();
    let propTime = new Date();
    /*
     * NOTE: camera formula sometimes is incorrect, but has been stable for over a year
     * NOTE: Looks wrong again as of 8/29/2020 - time of year issue?
     * NOTE: Could camera be related to daylight savings time? Subtracting one hour from selected date works
     */
    const doy = (0,_ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.getDayOfYear)(selectedDate);
    const modifier = 1000 * 60 * 60 * (-11.23 + 0.065666667 * doy);
    propTime.setUTCHours(selectedDate.getUTCHours());
    // + (selectedDate.getUTCMonth() * 2 - 11) / 2); // Offset has to account for time of year. Add 2 Hours per month into the year starting at -12.
    propTime.setUTCMinutes(selectedDate.getUTCMinutes());
    propTime.setUTCSeconds(selectedDate.getUTCSeconds());
    propTime = new Date(propTime.getTime() * 1 + modifier);
    realTime.setUTCHours(0, 0, 0, 0);
    const longOffset = (((propTime.getTime() - realTime.getTime()) / 60 / 60 / 1000) % 24) * 15; // 15 Degress Per Hour longitude Offset
    return normalizeAngle(((lon + longOffset) * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD));
};
/**
 * This function converts latitude in degrees to pitch in radians.
 */
const lat2pitch = (lat) => {
    const QUARTER_TAU = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU / 4;
    let pitch = lat * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD;
    pitch = Math.min(Math.max(pitch, -QUARTER_TAU), QUARTER_TAU);
    return pitch;
};
/**
 * This function converts altitude in kilometers to a zoom level between MIN_ZOOM_LEVEL and MAX_ZOOM_LEVEL for the camera.
 */
const alt2zoom = (alt, minZoomDistance, maxZoomDistance, minDistanceFromSatellite) => {
    if (minZoomDistance > maxZoomDistance) {
        throw new Error('minZoomDistance must be less than maxZoomDistance');
    }
    const distanceFromCenter = alt + _constants__WEBPACK_IMPORTED_MODULE_1__.RADIUS_OF_EARTH + minDistanceFromSatellite;
    const zoomLevel = ((distanceFromCenter - minZoomDistance) / (maxZoomDistance - minZoomDistance)) ** (1 / _constants__WEBPACK_IMPORTED_MODULE_1__.ZOOM_EXP);
    return Number.isNaN(zoomLevel) ? 0.5 : Math.min(Math.max(zoomLevel, 0), 1);
};
const isLeapYear_ = (dateIn) => {
    const year = dateIn.getUTCFullYear();
    if ((year & 3) !== 0) {
        return false;
    }
    return year % 100 !== 0 || year % 400 === 0;
};
/**
 * Calculate the julian day from a calendar date.
 */
const jday = (year, month, day, hour, minute, second) => {
    // Any negative values throw an error
    if (year < 0 || month < 1 || day < 0 || hour < 0 || minute < 0 || second < 0) {
        throw new Error('Invalid negative value');
    }
    // Validate month
    if (month > 12) {
        throw new Error('Invalid month value');
    }
    // Validate day
    if (day > 31) {
        throw new Error('Invalid day value');
    }
    // Validate hour
    if (hour > 23) {
        throw new Error('Invalid hour value');
    }
    // Validate minute
    if (minute > 59) {
        throw new Error('Invalid minute value');
    }
    // Validate second
    if (second > 60) {
        throw new Error('Invalid second value');
    }
    return (367.0 * year -
        Math.trunc(7 * (year + Math.trunc((month + 9) / 12.0)) * 0.25) +
        Math.trunc((275 * month) / 9.0) +
        day +
        1721013.5 +
        ((second / 60.0 + minute) / 60.0 + hour) / 24.0);
};
/**
 * Converts a local date to a zulu (UTC) date.
 */
const localToZulu = (date) => new Date(date.toISOString());
const dateFromJday = (year, day) => {
    if (year < 0) {
        throw new Error('Invalid negative value');
    }
    if (day < 1 || day > 366) {
        throw new Error('Invalid day value');
    }
    const date = new Date(Date.UTC(year, 0)); // initialize a date in `year-01-01` in UTC
    if (isLeapYear_(date)) {
        if (day > 366) {
            throw new Error('Invalid day value');
        }
    }
    else if (day > 365) {
        throw new Error('Invalid day value');
    }
    return new Date(date.setUTCDate(day)); // set the UTC date to the specified day
};
const dateToLocalInIso = (date) => {
    const offsetMs = -date.getTimezoneOffset() * 60 * 1000;
    const localDate = new Date(date.getTime() + offsetMs);
    const iso = localDate.toISOString().replace('T', ' ');
    return `${iso.slice(0, 19)} ${iso.slice(25, 31)}`;
};


}),
"./src/webworker/orbit-cruncher-interfaces.ts": 
/*!****************************************************!*\
  !*** ./src/webworker/orbit-cruncher-interfaces.ts ***!
  \****************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  OrbitCruncherMsgType: () => (OrbitCruncherMsgType),
  OrbitDrawTypes: () => (OrbitDrawTypes)
});
/**
 * Message Types for communication with the Orbit Cruncher Worker
 * These must match those in orbitCruncher.ts
 */
var OrbitCruncherMsgType;
(function (OrbitCruncherMsgType) {
    OrbitCruncherMsgType[OrbitCruncherMsgType["INIT"] = 0] = "INIT";
    OrbitCruncherMsgType[OrbitCruncherMsgType["UPDATE"] = 1] = "UPDATE";
    OrbitCruncherMsgType[OrbitCruncherMsgType["CHANGE_ORBIT_TYPE"] = 2] = "CHANGE_ORBIT_TYPE";
    OrbitCruncherMsgType[OrbitCruncherMsgType["MISSILE_UPDATE"] = 3] = "MISSILE_UPDATE";
    OrbitCruncherMsgType[OrbitCruncherMsgType["SATELLITE_UPDATE"] = 4] = "SATELLITE_UPDATE";
    OrbitCruncherMsgType[OrbitCruncherMsgType["SETTINGS_UPDATE"] = 5] = "SETTINGS_UPDATE";
    OrbitCruncherMsgType[OrbitCruncherMsgType["RESPONSE_READY"] = 6] = "RESPONSE_READY";
    OrbitCruncherMsgType[OrbitCruncherMsgType["RESPONSE_DATA"] = 7] = "RESPONSE_DATA";
})(OrbitCruncherMsgType || (OrbitCruncherMsgType = {}));
var OrbitDrawTypes;
(function (OrbitDrawTypes) {
    OrbitDrawTypes[OrbitDrawTypes["ORBIT"] = 0] = "ORBIT";
    OrbitDrawTypes[OrbitDrawTypes["TRAIL"] = 1] = "TRAIL";
})(OrbitDrawTypes || (OrbitDrawTypes = {}));


}),
"./src/webworker/positionCruncher/calculations.ts": 
/*!********************************************************!*\
  !*** ./src/webworker/positionCruncher/calculations.ts ***!
  \********************************************************/
(function (__unused_webpack_module, __webpack_exports__, __webpack_require__) {
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  checkSunExclusion: () => (checkSunExclusion),
  createLatLonAlt: () => (createLatLonAlt),
  createLatLonAltRad: () => (createLatLonAltRad),
  isInFov: () => (isInFov),
  isInValidElevation: () => (isInValidElevation),
  propTime: () => (propTime),
  setupTimeVariables: () => (setupTimeVariables)
});
/* ESM import */var _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @ootk/src/main */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../../engine/utils/external/meuusjs */ "./src/engine/utils/external/meuusjs.ts");
/* ESM import */var _engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../../engine/utils/transforms */ "./src/engine/utils/transforms.ts");



/* Returns Current Propagation Time */
const propTime = (dynamicOffsetEpoch, staticOffset, propRate) => {
    const now = new Date();
    const dynamicPropOffset = now.getTime() - dynamicOffsetEpoch;
    now.setTime(dynamicOffsetEpoch + staticOffset + dynamicPropOffset * propRate);
    return now;
};
const checkSunExclusion = (sensor, j, gmst, now) => {
    const jdo = new _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.JulianDay(j); // now
    // eslint-disable-next-line new-cap
    const coord = _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.EclCoordfromWgs84(0, 0, 0);
    // eslint-disable-next-line new-cap
    const coord2 = _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.EclCoordfromWgs84(sensor.lat, sensor.lon, sensor.alt);
    // AZ / EL Calculation
    const tp = _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.Solar.topocentricPosition(jdo, coord, false);
    const tpRel = _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.Solar.topocentricPosition(jdo, coord2, false);
    const sunAz = (tp.hz.az * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG + (180 % 360));
    const sunEl = ((tp.hz.alt * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG) % 360);
    const sunElRel = (tpRel.hz.alt * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG) % 360;
    // Range Calculation
    const T = new _engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.JulianDay(_engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.JulianDay.dateToJD(now)).jdJ2000Century();
    let sunG = (_engine_utils_external_meuusjs__WEBPACK_IMPORTED_MODULE_1__.A.Solar.meanAnomaly(T) * 180) / _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.PI;
    sunG %= 360.0;
    const sunR = 1.00014 - 0.01671 * Math.cos(sunG) - 0.00014 * Math.cos(2 * sunG);
    const sunRange = ((sunR * 149597870700) / 1000); // au to km conversion
    // RAE to ECI
    const sunECI = (0,_ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.rae2eci)({ az: sunAz, el: sunEl, rng: sunRange }, { lat: 0, lon: 0, alt: 0 }, gmst);
    return sensor && (sensor.type === _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.OPTICAL || sensor.type === _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.SpaceObjectType.OBSERVER) && sunElRel > -6 ? [true, sunECI] : [false, sunECI];
};
const isInFov = (sensor, lookangles) => {
    if (!lookangles) {
        return 0;
    }
    const { az, el, rng } = lookangles;
    sensor.minAz2 ??= Infinity;
    sensor.maxAz2 ??= -Infinity;
    sensor.minEl2 ??= Infinity;
    sensor.maxEl2 ??= -Infinity;
    sensor.minRng2 ??= Infinity;
    sensor.maxRng2 ??= -Infinity;
    if (sensor.minAz > sensor.maxAz) {
        if (((az >= sensor.minAz || az <= sensor.maxAz) && el >= sensor.minEl && el <= sensor.maxEl && rng <= sensor.maxRng && rng >= sensor.minRng) ||
            ((az >= (sensor.minAz2) || az <= sensor.maxAz2) && el >= sensor.minEl2 && el <= sensor.maxEl2 && rng <= sensor.maxRng2 && rng >= sensor.minRng2)) {
            return 1;
        }
    }
    else if ((az >= sensor.minAz && az <= sensor.maxAz && el >= sensor.minEl && el <= sensor.maxEl && rng <= sensor.maxRng && rng >= sensor.minRng) ||
        (az >= sensor.minAz2 && az <= sensor.maxAz2 && el >= sensor.minEl2 && el <= sensor.maxEl2 && rng <= sensor.maxRng2 && rng >= sensor.minRng2)) {
        return 1;
    }
    return 0;
};
const setupTimeVariables = (dynamicOffsetEpoch, staticOffset, propRate, isSunlightView, sensors) => {
    const now = propTime(dynamicOffsetEpoch, staticOffset, propRate);
    const j = (0,_engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__.jday)(now.getUTCFullYear(), now.getUTCMonth() + 1, // Note, this function requires months in range 1-12.
    now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds()) +
        now.getUTCMilliseconds() * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.MILLISECONDS_TO_DAYS;
    const gmst = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.gstime(j);
    let isSunExclusion = false;
    let sunEci = { x: 0, y: 0, z: 0 };
    if (isSunlightView && sensors?.length === 1) {
        // TODO: Sun exclusion should be calculated for each sensor
        [isSunExclusion, sunEci] = checkSunExclusion(sensors[0], j, gmst, now);
    }
    const j2 = (0,_engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__.jday)(now.getUTCFullYear(), now.getUTCMonth() + 1, // Note, this function requires months in range 1-12.
    now.getUTCDate(), now.getUTCHours(), now.getUTCMinutes(), now.getUTCSeconds() + 1) +
        now.getUTCMilliseconds() * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.MILLISECONDS_TO_DAYS;
    const gmstNext = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.gstime(j2);
    return {
        now,
        j,
        gmst,
        gmstNext,
        isSunExclusion,
        sunEci,
    };
};
const createLatLonAltRad = (lat, lon, alt) => ({
    lon,
    lat,
    alt,
});
const createLatLonAlt = (lat, lon, alt) => ({
    lon: (lon * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
    lat: (lat * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.RAD2DEG),
    alt,
});
const isInValidElevation = (rae, selectedSatFOV) => rae.el > 90 - selectedSatFOV;


}),

});
/************************************************************************/
// The module cache
var __webpack_module_cache__ = {};

// The require function
function __webpack_require__(moduleId) {

// Check if module is in cache
var cachedModule = __webpack_module_cache__[moduleId];
if (cachedModule !== undefined) {
return cachedModule.exports;
}
// Create a new module (and put it into the cache)
var module = (__webpack_module_cache__[moduleId] = {
exports: {}
});
// Execute the module function
__webpack_modules__[moduleId](module, module.exports, __webpack_require__);

// Return the exports of the module
return module.exports;

}

/************************************************************************/
// webpack/runtime/define_property_getters
(() => {
__webpack_require__.d = (exports, definition) => {
	for(var key in definition) {
        if(__webpack_require__.o(definition, key) && !__webpack_require__.o(exports, key)) {
            Object.defineProperty(exports, key, { enumerable: true, get: definition[key] });
        }
    }
};
})();
// webpack/runtime/has_own_property
(() => {
__webpack_require__.o = (obj, prop) => (Object.prototype.hasOwnProperty.call(obj, prop))
})();
// webpack/runtime/make_namespace_object
(() => {
// define __esModule on exports
__webpack_require__.r = (exports) => {
	if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
		Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
	}
	Object.defineProperty(exports, '__esModule', { value: true });
};
})();
// webpack/runtime/rspack_version
(() => {
__webpack_require__.rv = () => ("1.4.11")
})();
// webpack/runtime/rspack_unique_id
(() => {
__webpack_require__.ruid = "bundler=rspack@1.4.11";

})();
/************************************************************************/
var __webpack_exports__ = {};
// This entry needs to be wrapped in an IIFE because it needs to be isolated against other modules in the chunk.
(() => {

/*!****************************************!*\
  !*** ./src/webworker/orbitCruncher.ts ***!
  \****************************************/
__webpack_require__.r(__webpack_exports__);
__webpack_require__.d(__webpack_exports__, {
  onMessage: () => (onMessage)
});
/* ESM import */var _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(/*! @ootk/src/main */ "./src/engine/ootk/src/main.ts");
/* ESM import */var _engine_utils_constants__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(/*! ../engine/utils/constants */ "./src/engine/utils/constants.ts");
/* ESM import */var _engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(/*! ../engine/utils/transforms */ "./src/engine/utils/transforms.ts");
/* ESM import */var _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(/*! ./orbit-cruncher-interfaces */ "./src/webworker/orbit-cruncher-interfaces.ts");
/* ESM import */var _positionCruncher_calculations__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(/*! ./positionCruncher/calculations */ "./src/webworker/positionCruncher/calculations.ts");





let dynamicOffsetEpoch;
let staticOffset = 0;
let propRate = 1.0;
/** CONSTANTS */
const objCache = [];
let numberOfSegments;
let orbitType = _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitDrawTypes.ORBIT;
let orbitFadeFactor = 1.0;
let numberOfOrbitsToDraw = 1;
const trapIfJest_ = (cb) => {
    try {
        cb();
    }
    catch (e) {
        // If Jest isn't running then throw the error
        if (!process) {
            throw e;
        }
    }
};
const onMessage = (m) => {
    switch (m.data.type) {
        case _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.INIT:
            handleMsgInit_(m.data);
            return;
        case _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.SATELLITE_UPDATE:
            handleMsgSatelliteUpdate_(m.data);
            break;
        case _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.MISSILE_UPDATE:
            handleMsgMissileUpdate_(m.data);
            break;
        case _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.SETTINGS_UPDATE:
            handleMsgSettingsUpdate_(m.data);
            break;
        case _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.CHANGE_ORBIT_TYPE:
            handleMsgChangeOrbitType(m.data);
            return;
        default:
            return;
    }
    if (m.data.type === _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.SATELLITE_UPDATE || m.data.type === _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.MISSILE_UPDATE) {
        updateOrbitData_(m.data);
    }
};
const updateOrbitData_ = (data) => {
    /*
    * TODO: figure out how to calculate the orbit points on constant
    * position slices, not timeslices (ugly perigees on HEOs)
    */
    dynamicOffsetEpoch = data.dynamicOffsetEpoch;
    staticOffset = data.staticOffset;
    propRate = data.propRate;
    const id = data.id;
    let isEcfOutput = data.isEcfOutput || false;
    const pointsOut = new Float32Array((numberOfSegments + 1) * 4);
    const len = numberOfSegments + 1;
    let i = 0;
    // Calculate Missile Orbits
    if (objCache[id].missile) {
        while (i < len) {
            const missile = objCache[id];
            if (missile.latList?.length === 0) {
                pointsOut[i * 4] = 0;
                pointsOut[i * 4 + 1] = 0;
                pointsOut[i * 4 + 2] = 0;
                pointsOut[i * 4 + 3] = 0;
                i++;
            }
            else {
                drawMissileSegment_(missile, i, pointsOut, len);
                i++;
            }
        }
    }
    else if (objCache[id].ignore || !objCache[id].satrec) {
        // Invalid objects or OemSatellite with no TLEs
        trapIfJest_(() => {
            postMessage({
                type: _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.RESPONSE_DATA,
                pointsOut,
                satId: id,
            });
        });
        return;
    }
    else {
        const nowDate = (0,_positionCruncher_calculations__WEBPACK_IMPORTED_MODULE_4__.propTime)(dynamicOffsetEpoch, staticOffset, propRate);
        const nowJ = (0,_engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__.jday)(nowDate.getUTCFullYear(), nowDate.getUTCMonth() + 1, nowDate.getUTCDate(), nowDate.getUTCHours(), nowDate.getUTCMinutes(), nowDate.getUTCSeconds()) +
            nowDate.getUTCMilliseconds() * 1.15741e-8; // days per millisecond
        const satelliteObject = objCache[id];
        const satrec = satelliteObject.satrec;
        const now = (nowJ - satrec.jdsatepoch) * 1440.0; // in minutes
        // Calculate Satellite Orbits
        const period = (2 * Math.PI) / satrec.no; // convert rads/min to min
        let timeslice = period / numberOfSegments;
        // If a ECF output and  Geostationary orbit, then we can draw multiple orbits
        if (isEcfOutput && period > 1420 && period < 1460 && satrec.ecco < 0.05) {
            timeslice *= numberOfOrbitsToDraw;
        }
        else {
            isEcfOutput = false;
        }
        if (orbitType === _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitDrawTypes.ORBIT) {
            while (i < len) {
                drawTleOrbitSegment_(now, i, timeslice, id, isEcfOutput, period, pointsOut, len);
                i++;
            }
        }
        else if (orbitType === _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitDrawTypes.TRAIL) {
            while (i < len) {
                drawTleOrbitSegmentTrail_(now, i, timeslice, id, isEcfOutput, period, pointsOut, len);
                i++;
            }
        }
    }
    // TODO: Explore SharedArrayBuffer Options
    trapIfJest_(() => {
        postMessage({
            type: _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.RESPONSE_DATA,
            pointsOut,
            satId: id,
        });
    });
};
const drawMissileSegment_ = (missile, i, pointsOut, len) => {
    const x = Math.round(missile.altList.length * (i / numberOfSegments));
    const missileTime = (0,_positionCruncher_calculations__WEBPACK_IMPORTED_MODULE_4__.propTime)(dynamicOffsetEpoch, staticOffset, propRate);
    const j = (0,_engine_utils_transforms__WEBPACK_IMPORTED_MODULE_2__.jday)(missileTime.getUTCFullYear(), missileTime.getUTCMonth() + 1, // Note, this function requires months in range 1-12.
    missileTime.getUTCDate(), missileTime.getUTCHours(), missileTime.getUTCMinutes(), missileTime.getUTCSeconds()) +
        missileTime.getUTCMilliseconds() * 1.15741e-8; // days per millisecond
    const gmst = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.gstime(j);
    const cosLat = Math.cos(missile.latList[x] * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    const sinLat = Math.sin(missile.latList[x] * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD);
    const cosLon = Math.cos(missile.lonList[x] * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD + gmst);
    const sinLon = Math.sin(missile.lonList[x] * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.DEG2RAD + gmst);
    pointsOut[i * 4] = (_engine_utils_constants__WEBPACK_IMPORTED_MODULE_1__.RADIUS_OF_EARTH + missile.altList[x]) * cosLat * cosLon;
    pointsOut[i * 4 + 1] = (_engine_utils_constants__WEBPACK_IMPORTED_MODULE_1__.RADIUS_OF_EARTH + missile.altList[x]) * cosLat * sinLon;
    pointsOut[i * 4 + 2] = (_engine_utils_constants__WEBPACK_IMPORTED_MODULE_1__.RADIUS_OF_EARTH + missile.altList[x]) * sinLat;
    pointsOut[i * 4 + 3] = Math.min(orbitFadeFactor * (len / (i + 1)), 1.0);
};
const drawTleOrbitSegmentTrail_ = (now, i, timeslice, id, isEcfOutput, period, pointsOut, len) => {
    const t = now + i * timeslice;
    const sv = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.propagate(objCache[id].satrec, t);
    if (!sv) {
        pointsOut[i * 4] = 0;
        pointsOut[i * 4 + 1] = 0;
        pointsOut[i * 4 + 2] = 0;
        pointsOut[i * 4 + 3] = 0;
        return;
    }
    let pos = sv.position;
    if (isEcfOutput) {
        pos = (0,_ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.eci2ecf)(pos, (i * timeslice * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU) / period);
    }
    pointsOut[i * 4] = pos.x;
    pointsOut[i * 4 + 1] = pos.y;
    pointsOut[i * 4 + 2] = pos.z;
    pointsOut[i * 4 + 3] = i < len / 40 ? Math.min(orbitFadeFactor * (len / 40 / (2 * (i + 1))), 1.0) : 0.0;
};
const drawTleOrbitSegment_ = (now, i, timeslice, id, isEcfOutput, period, pointsOut, len) => {
    const t = now + i * timeslice;
    const sv = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.propagate(objCache[id].satrec, t);
    if (!sv) {
        pointsOut[i * 4] = 0;
        pointsOut[i * 4 + 1] = 0;
        pointsOut[i * 4 + 2] = 0;
        pointsOut[i * 4 + 3] = 0;
        return;
    }
    let pos = sv.position;
    if (isEcfOutput) {
        pos = (0,_ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.eci2ecf)(pos, (i * timeslice * _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.TAU) / period);
    }
    pointsOut[i * 4] = pos.x;
    pointsOut[i * 4 + 1] = pos.y;
    pointsOut[i * 4 + 2] = pos.z;
    pointsOut[i * 4 + 3] = Math.min(orbitFadeFactor * (len / (i + 1)), 1.0);
};
const handleMsgInit_ = (data) => {
    orbitFadeFactor = data.orbitFadeFactor ?? 1.0;
    numberOfOrbitsToDraw = data.numberOfOrbitsToDraw ?? 1;
    numberOfSegments = data.numSegs;
    const objData = JSON.parse(data.objData);
    const sLen = objData.length - 1;
    let i = -1;
    while (i < sLen) {
        i++;
        if (objData[i].missile) {
            objCache[i] = objData[i];
        }
        else if (objData[i].ignore) {
            objCache[i] = { ignore: true };
        }
        else if (objData[i].tle1 && objData[i].tle2) {
            objCache[i] = {
                satrec: _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.createSatrec(objData[i].tle1, objData[i].tle2),
            };
        }
        else {
            throw new Error('Invalid Object Data');
        }
    }
    postMessage({
        type: _orbit_cruncher_interfaces__WEBPACK_IMPORTED_MODULE_3__.OrbitCruncherMsgType.RESPONSE_READY,
    });
};
const handleMsgSatelliteUpdate_ = (data) => {
    // If new orbit
    if (data.tle1 && data.tle2) {
        const satelliteCacheEntry = objCache[data.id];
        satelliteCacheEntry.satrec = _ootk_src_main__WEBPACK_IMPORTED_MODULE_0__.Sgp4.createSatrec(data.tle1, data.tle2);
    }
};
const handleMsgMissileUpdate_ = (data) => {
    if (data.latList && data.lonList && data.altList) {
        const missileCacheEntry = objCache[data.id];
        missileCacheEntry.latList = data.latList;
        missileCacheEntry.lonList = data.lonList;
        missileCacheEntry.altList = data.altList;
    }
};
const handleMsgSettingsUpdate_ = (data) => {
    numberOfOrbitsToDraw = data.numberOfOrbitsToDraw ?? numberOfOrbitsToDraw;
};
const handleMsgChangeOrbitType = (data) => {
    orbitType = data.orbitType;
};
// Set up the web worker
onmessage = (m) => {
    trapIfJest_(() => onMessage(m));
};

})();

})()
;
//# sourceMappingURL=orbitCruncher.js.map