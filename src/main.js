import * as THREE from "three";
import { OrbitControls } from 'jsm/controls/OrbitControls.js';

import getStarfield from "./getStarfield.js";
import { getFresnelMat } from "./getFresnelMat.js";

// Debug: Check if Three.js loaded
console.log("Three.js loaded:", THREE);
console.log("OrbitControls loaded:", OrbitControls);

/* ---------------------- Scene & Earth (kept from your original) ---------------------- */
const w = window.innerWidth;
const h = window.innerHeight;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, w / h, 0.1, 1000);
camera.position.set(0, 2.5, 8);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(w, h);
document.body.appendChild(renderer.domElement);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
renderer.toneMappingExposure = 1.2;

const earthGroup = new THREE.Group();
earthGroup.rotation.z = -23.4 * Math.PI / 180;
scene.add(earthGroup);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.maxDistance = 250;
controls.target.set(5, 0, 0);
controls.update();

const detail = 12;
const loader = new THREE.TextureLoader();

// Add error handling for texture loading
loader.load("./textures/earthmap10k.jpg", 
  (texture) => console.log("Earth map loaded successfully"),
  undefined,
  (error) => console.error("Error loading earth map:", error)
);

const geometry = new THREE.IcosahedronGeometry(1, detail);
const material = new THREE.MeshPhongMaterial({
  map: loader.load("./textures/earthmap10k.jpg"),
  specularMap: loader.load("./textures/earthspec10k.jpg"),
  bumpMap: loader.load("/textures/earthbump10k.jpg"),
  bumpScale: 0.04,
});
const earthMesh = new THREE.Mesh(geometry, material);
earthGroup.add(earthMesh);

const lightsMat = new THREE.MeshBasicMaterial({
  map: loader.load("./textures/earthlights10k.jpg"),
  blending: THREE.AdditiveBlending,
});
const lightsMesh = new THREE.Mesh(geometry, lightsMat);
earthGroup.add(lightsMesh);

const cloudsMat = new THREE.MeshStandardMaterial({
  map: loader.load("./textures/earthcloudmap.jpg"),
  transparent: true,
  opacity: 0.4,
  blending: THREE.AdditiveBlending,
  alphaMap: loader.load('./textures/earthcloudmaptrans.jpg'),
});
const cloudsMesh = new THREE.Mesh(geometry, cloudsMat);
cloudsMesh.scale.setScalar(1.003);
earthGroup.add(cloudsMesh);

const fresnelMat = getFresnelMat();
const glowMesh = new THREE.Mesh(geometry, fresnelMat);
glowMesh.scale.setScalar(1.01);
earthGroup.add(glowMesh);

const stars = getStarfield({numStars: 2000});
scene.add(stars);

const sunLight = new THREE.DirectionalLight(0xffffff, 2.2);
sunLight.position.set(-2, 0.5, 1.5);
scene.add(sunLight);
// Add gentle ambient to brighten overall scene
scene.add(new THREE.AmbientLight(0xffffff, 0.6));
// Ensure initial view frames Earth
camera.lookAt(5, 0, 0);

/* ---------------------- Simple Sun + Earth orbit positioning ---------------------- */
/* We'll place a small sun at origin. We'll render Earth's orbit ring at 1 AU (scaled). */
const AU_TO_SCENE = 5; // 1 AU == 5 scene units (tweak as you like)
// Create a more realistic sun with emissive material
const sunMat = new THREE.MeshBasicMaterial({ 
  color: 0xffaa44,
  emissive: 0xff6600,
  emissiveIntensity: 0.8
});
const sunMesh = new THREE.Mesh(new THREE.SphereGeometry(0.5, 32, 32), sunMat);
scene.add(sunMesh);

// Add a subtle glow using lensflare textures
const glowTex = new THREE.TextureLoader().load('./textures/lensflare0.png');
const sunGlow = new THREE.Sprite(new THREE.SpriteMaterial({ 
  map: glowTex, 
  color: 0xffcc88, 
  transparent: true, 
  opacity: 0.6, 
  depthWrite: false,
  blending: THREE.AdditiveBlending
}));
sunGlow.scale.set(2.2, 2.2, 1);
sunMesh.add(sunGlow);

// Add additional lensflare effect for more depth
const lensflareTex = new THREE.TextureLoader().load('./textures/lensflare3.png');
const sunLensflare = new THREE.Sprite(new THREE.SpriteMaterial({ 
  map: lensflareTex, 
  color: 0xffaa44, 
  transparent: true, 
  opacity: 0.4, 
  depthWrite: false,
  blending: THREE.AdditiveBlending
}));
sunLensflare.scale.set(3.0, 3.0, 1);
sunMesh.add(sunLensflare);

/* Move the Earth mesh to its orbital position at 1 AU (so orbits appear relative to sun) */
const earthOrbitDistance = 1 * AU_TO_SCENE;
earthGroup.position.set(earthOrbitDistance, 0, 0);
earthMesh.scale.setScalar(0.25); // make planet smaller relative to orbit visualization
lightsMesh.scale.copy(earthMesh.scale);
cloudsMesh.scale.setScalar(1.003 * 0.25);
glowMesh.scale.setScalar(1.01 * 0.25);
// Make Earth bigger and glow brighter
earthMesh.scale.setScalar(0.4);
lightsMesh.scale.copy(earthMesh.scale);
cloudsMesh.scale.setScalar(1.003 * 0.4);
glowMesh.scale.setScalar(1.01 * 0.4);
glowMesh.material.uniforms.uOpacity.value = 0.35;

/* Draw Earth's orbit ring (for reference) */
function drawOrbitRing(radius) {
  const segments = 360;
  const points = [];
  for (let i = 0; i <= segments; i++) {
    const theta = (i / segments) * Math.PI * 2;
    points.push(new THREE.Vector3(Math.cos(theta) * radius, 0, Math.sin(theta) * radius));
  }
  const g = new THREE.BufferGeometry().setFromPoints(points);
  const mat = new THREE.LineBasicMaterial({ color: 0x4444ff, transparent: true, opacity: 0.25 });
  const l = new THREE.LineLoop(g, mat);
  scene.add(l);
}
drawOrbitRing(earthOrbitDistance);

/* ---------------------- Orbit visualization helpers (Kepler -> Cartesian) ---------------------- */

/* Utility: solve Kepler's equation M = E - e*sin E for E given M & e */
function solveKepler(M, e, maxIter = 60, tol = 1e-6) {
  // normalize M to -PI..PI
  M = ((M + Math.PI) % (2 * Math.PI)) - Math.PI;
  let E = e < 0.8 ? M : Math.PI;
  for (let i = 0; i < maxIter; i++) {
    const f = E - e * Math.sin(E) - M;
    const fp = 1 - e * Math.cos(E);
    const dE = f / fp;
    E -= dE;
    if (Math.abs(dE) < tol) break;
  }
  return E;
}

/* Create orbit points (in scene units) from classical orbital elements (all angles in radians) */
/* a in AU, e, i (rad), Omega (ra), omega (arg per), samples */
function orbitPointsFromElements(a, e, i, Omega, omega, samples = 360) {
  const points = [];
  for (let si = 0; si <= samples; si++) {
    const M = (si / samples) * 2 * Math.PI;
    const E = solveKepler(M, e);
    const nu = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
    const r = a * (1 - e * Math.cos(E)); // in AU
    // position in orbital plane
    const xOrb = r * Math.cos(nu);
    const yOrb = r * Math.sin(nu);
    // rotate to ecliptic / inertial frame via: r_vec = Rz(Omega) * Rx(i) * Rz(omega) * [xOrb, yOrb, 0]
    const cosO = Math.cos(Omega), sinO = Math.sin(Omega);
    const cosi = Math.cos(i), sini = Math.sin(i);
    const cosw = Math.cos(omega), sinw = Math.sin(omega);

    // position vector in ecliptic coordinates
    const X = (cosO * cosw - sinO * sinw * cosi) * xOrb + (-cosO * sinw - sinO * cosw * cosi) * yOrb;
    const Y = (sinO * cosw + cosO * sinw * cosi) * xOrb + (-sinO * sinw + cosO * cosw * cosi) * yOrb;
    const Z = (sinw * sini) * xOrb + (cosw * sini) * yOrb;

    points.push(new THREE.Vector3(X * AU_TO_SCENE, Z * AU_TO_SCENE, Y * AU_TO_SCENE));
  }
  return points;
}

/* Support multiple asteroids and projectiles */
const asteroidObjects = new Map(); // id -> { line, marker, points, t, color, meta }
const projectiles = new Set(); // active dropping asteroids
let lastSelectedAsteroidId = null;

function createAsteroidVisual(points, color = 0xffaa33) {
  const g = new THREE.BufferGeometry().setFromPoints(points);
  const line = new THREE.Line(g, new THREE.LineBasicMaterial({ color, linewidth: 2 }));
  const marker = new THREE.Mesh(
    new THREE.SphereGeometry(0.12, 12, 12),
    new THREE.MeshStandardMaterial({ color: 0xffaa66, emissive: 0x552200, emissiveIntensity: 1.5 })
  );
  scene.add(line);
  scene.add(marker);
  return { line, marker };
}

function randomNiceColor(seed = Math.random()) {
  const hue = (seed * 360) % 360;
  const color = new THREE.Color().setHSL(hue / 360, 0.7, 0.6);
  return color.getHex();
}

// Expose hooks for index.html (multi-add/remove/clear)
window.addAsteroidOrbit = function(orbitalData, id = undefined, meta = undefined) {
  if (!orbitalData) return;
  const key = id || `${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;
  if (asteroidObjects.has(key)) return; // already added
  const a = parseFloat(orbitalData.semi_major_axis) || parseFloat(orbitalData.semi_major_axis_au) || 1.5; // AU
  const e = parseFloat(orbitalData.eccentricity) || 0.2;
  const i = THREE.MathUtils.degToRad(parseFloat(orbitalData.inclination) || 0);
  const Omega = THREE.MathUtils.degToRad(parseFloat(orbitalData.ascending_node_longitude) || 0);
  const omega = THREE.MathUtils.degToRad(parseFloat(orbitalData.perihelion_argument) || 0);
  const pts = orbitPointsFromElements(a, e, i, Omega, omega, 1400);
  const color = randomNiceColor();
  const { line, marker } = createAsteroidVisual(pts, color);
  marker.userData.id = key;
  marker.userData.meta = meta;
  asteroidObjects.set(key, { line, marker, points: pts, t: Math.random(), color, meta });
  lastSelectedAsteroidId = key;
  return key;
};

window.removeAsteroid = function(id) {
  const obj = asteroidObjects.get(id);
  if (!obj) return;
  scene.remove(obj.line); obj.line.geometry.dispose(); obj.line.material.dispose();
  scene.remove(obj.marker); obj.marker.geometry.dispose(); obj.marker.material.dispose();
  asteroidObjects.delete(id);
};

window.clearAsteroids = function() {
  for (const [id, obj] of asteroidObjects.entries()) {
    scene.remove(obj.line); obj.line.geometry.dispose(); obj.line.material.dispose();
    scene.remove(obj.marker); obj.marker.geometry.dispose(); obj.marker.material.dispose();
  }
  asteroidObjects.clear();
};

// Click picking to show asteroid info panel when clicking a marker
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

function onClick(event) {
  const rect = renderer.domElement.getBoundingClientRect();
  mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  raycaster.setFromCamera(mouse, camera);
  const markers = Array.from(asteroidObjects.values()).map(o => o.marker);
  const hits = raycaster.intersectObjects(markers, false);
  if (hits.length) {
    const hit = hits[0].object;
    const meta = hit.userData.meta;
    const id = hit.userData.id;
    lastSelectedAsteroidId = id;
    const panel = document.getElementById("infoPanel");
    if (panel) {
      if (meta) {
        panel.innerHTML = `
          <h3>${meta.name}</h3>
          <p><b>ID:</b> ${meta.id}</p>
          <p><b>Size:</b> ${meta.est_diameter || 'N/A'} m</p>
          <p><b>Speed:</b> ${meta.speed || 'N/A'} km/h</p>
          <p><b>Nearest to Earth:</b> ${meta.nearest_distance_km || 'N/A'} km</p>
          <p><b>Discovered:</b> ${meta.discovery_date || 'N/A'}</p>
        `;
      } else {
        panel.innerHTML = `
          <h3>Asteroid ${id}</h3>
          <p>No additional metadata available.</p>
        `;
      }
      panel.style.display = "block";
    }
  }
}
renderer.domElement.addEventListener('click', onClick);

// Convert city lat/lon to world point for targeting (duplicate for local use)
function latLonToPoint(latDeg, lonDeg, radius = 0.4) {
  const lat = THREE.MathUtils.degToRad(latDeg);
  const lon = THREE.MathUtils.degToRad(lonDeg);
  const x = Math.cos(lat) * Math.cos(lon) * radius;
  const y = Math.sin(lat) * radius;
  const z = Math.cos(lat) * Math.sin(lon) * radius;
  return new THREE.Vector3(x, y, z).add(earthGroup.position);
}

// Drop straight-line projectile toward city; no orbit line
window.dropAsteroidTo = function(lat, lon, diameterMeters, velocityKmS, cityName) {
  window.clearAsteroids();
  const start = new THREE.Vector3(-10, 5, -10);
  const target = latLonToPoint(lat, lon, 0.4);
  const dir = target.clone().sub(start).normalize();
  const speed = Math.max(0.5, velocityKmS / 20);
  const size = THREE.MathUtils.clamp(diameterMeters / 1000, 0.06, 0.6);
  const marker = new THREE.Mesh(
    new THREE.SphereGeometry(size, 14, 14),
    new THREE.MeshStandardMaterial({ color: 0xff9966, emissive: 0x552200, emissiveIntensity: 2.0 })
  );
  marker.position.copy(start);
  scene.add(marker);
  const trail = new THREE.Line(
    new THREE.BufferGeometry().setFromPoints([start.clone(), start.clone()]),
    new THREE.LineBasicMaterial({ color: 0xff6633, transparent: true, opacity: 0.6 })
  );
  scene.add(trail);
  const projectile = { marker, trail, dir, speed, target, cityName };
  projectiles.add(projectile);
};
// Enhanced mitigation effects
window.applyMitigation = function(type) {
  // operate on last selected asteroid
  const ids = Array.from(asteroidObjects.keys());
  const lastId = ids[ids.length - 1];
  const obj = asteroidObjects.get(lastId);
  if (!obj || !obj.meta || !obj.meta.orbital_data) return;
  
  const od = obj.meta.orbital_data;
  const a0 = parseFloat(od.semi_major_axis) || parseFloat(od.semi_major_axis_au) || 1.5;
  const e0 = parseFloat(od.eccentricity) || 0.2;
  const i0 = THREE.MathUtils.degToRad(parseFloat(od.inclination) || 0);
  const Omega0 = THREE.MathUtils.degToRad(parseFloat(od.ascending_node_longitude) || 0);
  const omega0 = THREE.MathUtils.degToRad(parseFloat(od.perihelion_argument) || 0);
  
  let a = a0, e = e0, i = i0, Omega = Omega0, omega = omega0;
  let effectMessage = '';
  
  if (type === 'tractor') {
    // Gravity tractor: gradual orbital modification
    // Reduces eccentricity and increases semi-major axis (moves away from Earth)
    e = Math.max(0.01, e0 - 0.05); // Reduce eccentricity
    a = a0 + 0.1; // Increase semi-major axis
    // Slight inclination change
    i = i0 + THREE.MathUtils.degToRad(2);
    effectMessage = `Gravity Tractor Applied!\n\nOrbit modified:\n• Eccentricity: ${e0.toFixed(3)} → ${e.toFixed(3)}\n• Semi-major axis: ${a0.toFixed(3)} → ${a.toFixed(3)} AU\n• Inclination: ${(i0 * 180/Math.PI).toFixed(1)}° → ${(i * 180/Math.PI).toFixed(1)}°\n\nAsteroid trajectory deflected away from Earth.`;
  } else if (type === 'laser') {
    // Laser ablation: more aggressive orbital change
    // Changes perihelion argument and reduces semi-major axis
    omega = omega0 + THREE.MathUtils.degToRad(15); // Change argument of perihelion
    a = Math.max(0.8, a0 - 0.15); // Reduce semi-major axis
    e = Math.max(0.01, e0 + 0.03); // Slightly increase eccentricity
    // Significant inclination change
    i = i0 + THREE.MathUtils.degToRad(8);
    effectMessage = `Laser Ablation Applied!\n\nOrbit modified:\n• Perihelion argument: ${(omega0 * 180/Math.PI).toFixed(1)}° → ${(omega * 180/Math.PI).toFixed(1)}°\n• Semi-major axis: ${a0.toFixed(3)} → ${a.toFixed(3)} AU\n• Eccentricity: ${e0.toFixed(3)} → ${e.toFixed(3)}\n• Inclination: ${(i0 * 180/Math.PI).toFixed(1)}° → ${(i * 180/Math.PI).toFixed(1)}°\n\nAsteroid trajectory significantly altered.`;
  }
  
  // Generate new orbit points
  const pts = orbitPointsFromElements(a, e, i, Omega, omega, 1400);
  
  // Update visuals with color change to show modification
  const newColor = type === 'tractor' ? 0x00ff88 : 0xff8800; // Green for tractor, Orange for laser
  scene.remove(obj.line); 
  obj.line.geometry.dispose(); 
  obj.line.material.dispose();
  obj.line = new THREE.Line(
    new THREE.BufferGeometry().setFromPoints(pts), 
    new THREE.LineBasicMaterial({ color: newColor, linewidth: 3 })
  );
  scene.add(obj.line);
  obj.points = pts;
  obj.t = 0;
  
  // Update marker color
  obj.marker.material.color.setHex(newColor);
  
  // Write back orbital data
  od.semi_major_axis = a; 
  od.eccentricity = e;
  od.inclination = i * 180/Math.PI;
  od.ascending_node_longitude = Omega * 180/Math.PI;
  od.perihelion_argument = omega * 180/Math.PI;
  
  // Show effect message
  showImpactModal(effectMessage);
};

// Nudge currently selected asteroid by tweaking argument of perihelion a bit
window.nudgeSelectedAsteroid = function() {
  if (!lastSelectedAsteroidId) return;
  const obj = asteroidObjects.get(lastSelectedAsteroidId);
  if (!obj || !obj.meta || !obj.meta.orbital_data) return;
  const od = obj.meta.orbital_data;
  const a = parseFloat(od.semi_major_axis) || parseFloat(od.semi_major_axis_au) || 1.5;
  const e = parseFloat(od.eccentricity) || 0.2;
  const i = THREE.MathUtils.degToRad(parseFloat(od.inclination) || 0);
  const Omega = THREE.MathUtils.degToRad(parseFloat(od.ascending_node_longitude) || 0);
  const omegaDeg = (parseFloat(od.perihelion_argument) || 0) + 3; // apply +3°
  od.perihelion_argument = omegaDeg; // persist back to meta
  const omega = THREE.MathUtils.degToRad(omegaDeg);
  const pts = orbitPointsFromElements(a, e, i, Omega, omega, 1400);
  // Replace visuals
  scene.remove(obj.line); obj.line.geometry.dispose(); obj.line.material.dispose();
  obj.line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), new THREE.LineBasicMaterial({ color: obj.color, linewidth: 2 }));
  scene.add(obj.line);
  obj.points = pts;
  obj.t = 0;
};

/* ---------------------- Impact Effects System ---------------------- */
const impactEffects = new Set();
let cameraShakeIntensity = 0;
let cameraShakeDuration = 0;
const originalCameraPosition = camera.position.clone();

// Audio context for impact sounds (optional)
let audioContext = null;
let impactSoundBuffer = null;

function initAudio() {
  try {
    audioContext = new (window.AudioContext || window.webkitAudioContext)();
    createImpactSound();
  } catch (e) {
    console.log('Audio not supported or disabled');
  }
}

function createImpactSound() {
  if (!audioContext) return;
  
  const sampleRate = audioContext.sampleRate;
  const duration = 0.5; // seconds
  const buffer = audioContext.createBuffer(1, sampleRate * duration, sampleRate);
  const data = buffer.getChannelData(0);
  
  // Create a low-frequency rumble with some high-frequency crack
  for (let i = 0; i < data.length; i++) {
    const t = i / sampleRate;
    const rumble = Math.sin(2 * Math.PI * 60 * t) * Math.exp(-t * 3) * 0.3;
    const crack = (Math.random() - 0.5) * Math.exp(-t * 8) * 0.1;
    data[i] = rumble + crack;
  }
  
  impactSoundBuffer = buffer;
}

function playImpactSound(intensity = 1.0) {
  if (!audioContext || !impactSoundBuffer) return;
  
  try {
    const source = audioContext.createBufferSource();
    const gainNode = audioContext.createGain();
    
    source.buffer = impactSoundBuffer;
    gainNode.gain.value = Math.min(1.0, intensity * 0.5);
    
    source.connect(gainNode);
    gainNode.connect(audioContext.destination);
    
    source.start();
  } catch (e) {
    console.log('Could not play impact sound:', e);
  }
}

function createExplosionEffect(position, scale = 1.0) {
  const explosionGroup = new THREE.Group();
  explosionGroup.position.copy(position);
  scene.add(explosionGroup);

  // Main explosion sphere
  const explosionGeometry = new THREE.SphereGeometry(0.1, 16, 16);
  const explosionMaterial = new THREE.MeshBasicMaterial({
    color: 0xff6600,
    transparent: true,
    opacity: 0.8,
    blending: THREE.AdditiveBlending
  });
  const explosionMesh = new THREE.Mesh(explosionGeometry, explosionMaterial);
  explosionGroup.add(explosionMesh);

  // Shockwave ring
  const shockwaveGeometry = new THREE.RingGeometry(0.05, 0.15, 32);
  const shockwaveMaterial = new THREE.MeshBasicMaterial({
    color: 0xffff00,
    transparent: true,
    opacity: 0.6,
    side: THREE.DoubleSide,
    blending: THREE.AdditiveBlending
  });
  const shockwaveMesh = new THREE.Mesh(shockwaveGeometry, shockwaveMaterial);
  shockwaveMesh.rotation.x = Math.PI / 2;
  explosionGroup.add(shockwaveMesh);

  // Debris particles
  const debrisCount = 20;
  const debrisGeometry = new THREE.BufferGeometry();
  const debrisPositions = new Float32Array(debrisCount * 3);
  const debrisVelocities = [];
  
  for (let i = 0; i < debrisCount; i++) {
    const i3 = i * 3;
    debrisPositions[i3] = 0;
    debrisPositions[i3 + 1] = 0;
    debrisPositions[i3 + 2] = 0;
    
    debrisVelocities.push({
      x: (Math.random() - 0.5) * 0.3,
      y: Math.random() * 0.4 + 0.1,
      z: (Math.random() - 0.5) * 0.3
    });
  }
  
  debrisGeometry.setAttribute('position', new THREE.BufferAttribute(debrisPositions, 3));
  const debrisMaterial = new THREE.PointsMaterial({
    color: 0xff4400,
    size: 0.02,
    transparent: true,
    opacity: 0.8,
    blending: THREE.AdditiveBlending
  });
  const debrisPoints = new THREE.Points(debrisGeometry, debrisMaterial);
  explosionGroup.add(debrisPoints);

  const effect = {
    group: explosionGroup,
    explosion: explosionMesh,
    shockwave: shockwaveMesh,
    debris: debrisPoints,
    debrisVelocities: debrisVelocities,
    age: 0,
    maxAge: 3.0,
    scale: scale
  };

  impactEffects.add(effect);
  return effect;
}

function createEarthDamage(position, craterSize) {
  // Create a dark crater mark on Earth surface
  const craterGeometry = new THREE.CircleGeometry(craterSize, 32);
  const craterMaterial = new THREE.MeshBasicMaterial({
    color: 0x2a1810,
    transparent: true,
    opacity: 0.8,
    side: THREE.DoubleSide
  });
  const craterMesh = new THREE.Mesh(craterGeometry, craterMaterial);
  
  // Position crater on Earth surface
  const earthCenter = earthGroup.position;
  const direction = position.clone().sub(earthCenter).normalize();
  craterMesh.position.copy(earthCenter).add(direction.multiplyScalar(0.4));
  craterMesh.lookAt(earthCenter);
  craterMesh.rotation.z += Math.PI / 2;
  
  scene.add(craterMesh);

  // Add fire/smoke effect
  const fireGeometry = new THREE.ConeGeometry(craterSize * 0.8, craterSize * 0.5, 8);
  const fireMaterial = new THREE.MeshBasicMaterial({
    color: 0xff4400,
    transparent: true,
    opacity: 0.6,
    blending: THREE.AdditiveBlending
  });
  const fireMesh = new THREE.Mesh(fireGeometry, fireMaterial);
  fireMesh.position.copy(craterMesh.position);
  fireMesh.lookAt(camera.position);
  scene.add(fireMesh);

  // Animate fire
  const fireEffect = {
    mesh: fireMesh,
    age: 0,
    maxAge: 5.0,
    originalScale: fireMesh.scale.clone()
  };

  impactEffects.add(fireEffect);
  
  // Remove crater after some time
  setTimeout(() => {
    scene.remove(craterMesh);
    craterMesh.geometry.dispose();
    craterMesh.material.dispose();
  }, 10000);

  return craterMesh;
}

function updateImpactEffects(dt) {
  for (const effect of Array.from(impactEffects)) {
    effect.age += dt;
    
    if (effect.age >= effect.maxAge) {
      // Clean up expired effects
      if (effect.group) {
        scene.remove(effect.group);
        effect.group.traverse((child) => {
          if (child.geometry) child.geometry.dispose();
          if (child.material) child.material.dispose();
        });
      }
      if (effect.mesh) {
        scene.remove(effect.mesh);
        effect.mesh.geometry.dispose();
        effect.mesh.material.dispose();
      }
      impactEffects.delete(effect);
      continue;
    }

    const progress = effect.age / effect.maxAge;
    
    if (effect.explosion) {
      // Animate explosion
      const explosionScale = (1 + progress * 3) * effect.scale;
      effect.explosion.scale.setScalar(explosionScale);
      effect.explosion.material.opacity = (1 - progress) * 0.8;
      
      // Animate shockwave
      const shockwaveScale = (1 + progress * 5) * effect.scale;
      effect.shockwave.scale.setScalar(shockwaveScale);
      effect.shockwave.material.opacity = (1 - progress) * 0.6;
      
      // Animate debris
      if (effect.debris && effect.debrisVelocities) {
        const positions = effect.debris.geometry.attributes.position.array;
        for (let i = 0; i < effect.debrisVelocities.length; i++) {
          const i3 = i * 3;
          const vel = effect.debrisVelocities[i];
          positions[i3] += vel.x * dt;
          positions[i3 + 1] += vel.y * dt;
          positions[i3 + 2] += vel.z * dt;
          vel.y -= 0.5 * dt; // gravity
        }
        effect.debris.geometry.attributes.position.needsUpdate = true;
        effect.debris.material.opacity = (1 - progress) * 0.8;
      }
    }
    
    if (effect.mesh && effect.originalScale) {
      // Animate fire effect
      const fireScale = (1 - progress * 0.5) * effect.scale;
      effect.mesh.scale.copy(effect.originalScale).multiplyScalar(fireScale);
      effect.mesh.material.opacity = (1 - progress) * 0.6;
    }
  }
}

function addCameraShake(intensity = 0.1, duration = 0.5) {
  cameraShakeIntensity = intensity;
  cameraShakeDuration = duration;
}

function updateCameraShake(dt) {
  if (cameraShakeDuration > 0) {
    cameraShakeDuration -= dt;
    const shakeAmount = cameraShakeIntensity * (cameraShakeDuration / 0.5);
    
    camera.position.x = originalCameraPosition.x + (Math.random() - 0.5) * shakeAmount;
    camera.position.y = originalCameraPosition.y + (Math.random() - 0.5) * shakeAmount;
    camera.position.z = originalCameraPosition.z + (Math.random() - 0.5) * shakeAmount;
    
    if (cameraShakeDuration <= 0) {
      camera.position.copy(originalCameraPosition);
    }
  }
}

/* ---------------------- Animation loop ---------------------- */
let lastTS = performance.now();
function animate(ts) {
  requestAnimationFrame(animate);
  const dt = (ts - lastTS) / 1000;
  lastTS = ts;

  // spin earth components (kept from original)
  earthMesh.rotation.y += 0.002;
  lightsMesh.rotation.y += 0.002;
  cloudsMesh.rotation.y += 0.0023;
  glowMesh.rotation.y += 0.002;

  // revolve Earth around the Sun
  const revSpeed = 0.05; // radians per second (tweak)
  const angle = ts * 0.001 * revSpeed;
  earthGroup.position.set(Math.cos(angle) * earthOrbitDistance, 0, Math.sin(angle) * earthOrbitDistance);

  // Update impact effects
  updateImpactEffects(dt);
  
  // Update camera shake
  updateCameraShake(dt);
  
  // Rotate sun slowly
  sunMesh.rotation.y += 0.001;
  
  // Add subtle rotation to lensflare sprites for more dynamic effect
  if (sunGlow) {
    sunGlow.rotation.z += 0.0005;
  }
  if (sunLensflare) {
    sunLensflare.rotation.z -= 0.0003;
  }

  // animate all asteroids along their orbits
  for (const obj of asteroidObjects.values()) {
    obj.t = (obj.t + dt * 0.03) % 1; // speed multiplier
    const idx = Math.floor(obj.t * (obj.points.length - 1));
    obj.marker.position.copy(obj.points[idx]);
  }

  // animate projectiles toward target
  for (const p of Array.from(projectiles)) {
    const step = p.dir.clone().multiplyScalar(dt * p.speed);
    p.marker.position.add(step);
    // extend trail
    const geom = p.trail.geometry;
    const oldAttr = geom.getAttribute('position');
    const old = Array.from(oldAttr.array);
    old.push(p.marker.position.x, p.marker.position.y, p.marker.position.z);
    const newAttr = new THREE.Float32BufferAttribute(old, 3);
    geom.setAttribute('position', newAttr);
    geom.setDrawRange(0, old.length / 3);
    geom.attributes.position.needsUpdate = true;

    // impact check (distance to Earth's center <= radius -> impact)
    const earthCenter = earthGroup.position;
    const earthRadius = 0.4; // scene units
    const toCenter = p.marker.position.clone().sub(earthCenter).length();
    const hitGround = toCenter <= earthRadius + 0.1;
    const reachedTarget = p.marker.position.distanceTo(p.target) < 0.2;
    if (hitGround || reachedTarget) {
      // Calculate impact parameters
      const densEl = document.getElementById('densInput');
      const velEl = document.getElementById('velInput');
      const diaEl = document.getElementById('diaInput');
      const density = densEl ? parseFloat(densEl.value || '3000') : 3000;
      const vKmS = velEl ? parseFloat(velEl.value || '20') : 20;
      const dia = diaEl ? parseFloat(diaEl.value || '150') : 150;
      
      // Calculate energy and crater size
      const r = dia / 2;
      const volume = (4/3) * Math.PI * Math.pow(r, 3);
      const mass = density * volume;
      const v = vKmS * 1000;
      const energyJ = 0.5 * mass * v * v;
      const kt = energyJ / 4.184e9;
      const mt = kt / 1000;
      const craterM = 20 * Math.cbrt(Math.max(mt, 0));
      let deaths = Math.round(Math.min(50_000_000, (mt * 20000)));
      
      // Create impact effects
      if (hitGround) {
        const impactPosition = p.marker.position.clone();
        const impactScale = Math.min(3.0, Math.max(0.5, mt / 10)); // Scale based on energy
        
        // Create explosion effect
        createExplosionEffect(impactPosition, impactScale);
        
        // Create Earth surface damage
        const craterSize = Math.min(0.3, Math.max(0.05, craterM / 1000)); // Scale crater to scene units
        createEarthDamage(impactPosition, craterSize);
        
        // Add camera shake based on impact energy
        const shakeIntensity = Math.min(0.3, Math.max(0.05, mt / 50));
        const shakeDuration = Math.min(2.0, Math.max(0.5, mt / 25));
        addCameraShake(shakeIntensity, shakeDuration);
        
        // Play impact sound
        const soundIntensity = Math.min(1.0, Math.max(0.1, mt / 20));
        playImpactSound(soundIntensity);
        
        // Calculate detailed impact information
        const cityInfo = window.cityData ? window.cityData[p.cityName] : { population: 1000000, coastal: false };
        const fireballRadius = 0.5 * Math.cbrt(mt) * 1000; // meters
        const airblastRadius = 1.5 * Math.cbrt(mt) * 1000; // meters
        const earthquakeMagnitude = 4.8 + 0.67 * Math.log10(mt);
        const tsunamiHeight = cityInfo.coastal ? Math.max(0, 0.1 * Math.cbrt(mt)) : 0;
        
        // More realistic death calculation
        const population = cityInfo.population || 1000000;
        let realisticDeaths = 0;
        if (airblastRadius < 5000) {
          realisticDeaths = Math.min(population, population * 0.95);
        } else if (airblastRadius < 10000) {
          realisticDeaths = Math.min(population, population * 0.7);
        } else if (airblastRadius < 20000) {
          realisticDeaths = Math.min(population, population * 0.3);
        } else if (airblastRadius < 50000) {
          realisticDeaths = Math.min(population, population * 0.1);
        }
        
        // Add tsunami casualties if coastal
        if (cityInfo.coastal && tsunamiHeight > 0) {
          realisticDeaths += Math.min(population * 0.1, tsunamiHeight * 1000);
        }
        
        // Show enhanced impact modal with detailed information
        const survivalRate = ((1 - realisticDeaths/population) * 100).toFixed(1);
        showImpactModal(`💥 IMPACT DETECTED! 💥

📍 Location: ${p.cityName || 'target'}
⚡ Energy: ${mt.toFixed(2)} Mt (${kt.toFixed(0)} kt)
🕳️ Crater: ${craterM.toFixed(0)} m diameter
🔥 Fireball: ${fireballRadius.toFixed(0)} m radius
💨 Airblast: ${airblastRadius.toFixed(0)} m radius
🌍 Earthquake: ${earthquakeMagnitude.toFixed(1)} Richter
${cityInfo.coastal ? `🌊 Tsunami: ${tsunamiHeight.toFixed(1)} m height` : 'No tsunami (inland)'}

👥 Population: ${population.toLocaleString()}
💀 Estimated Deaths: ${Math.round(realisticDeaths).toLocaleString()}
✅ Survival Rate: ${survivalRate}%

Impact effects are now visible in the simulation.`);
      } else {
        // Missed Earth - smaller explosion in space
        createExplosionEffect(p.marker.position, 0.5);
        addCameraShake(0.02, 0.2);
        playImpactSound(0.1); // Small sound for miss
        showImpactModal(`Asteroid missed Earth by ${(toCenter - earthRadius).toFixed(2)} units.\nNo impact damage.`);
      }
      
      // Clean up projectile
      scene.remove(p.marker); p.marker.geometry.dispose(); p.marker.material.dispose();
      scene.remove(p.trail); p.trail.geometry.dispose(); p.trail.material.dispose();
      projectiles.delete(p);
    }
  }

  controls.update();
  renderer.render(scene, camera);
}
animate(performance.now());

/* handle resize */
function handleWindowResize () {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}
window.addEventListener('resize', handleWindowResize, false);

/* Initialize audio on first user interaction */
document.addEventListener('click', () => {
  if (!audioContext) {
    initAudio();
  }
}, { once: true });

/* Clean up on navigation (optional) */
window.addEventListener("beforeunload", () => {
  renderer.dispose();
});

/* ---------------------- City focus and impact corridor ---------------------- */
const cityMarkers = new Map();
const corridorGroup = new THREE.Group();
scene.add(corridorGroup);

function latLonToVector3(latDeg, lonDeg, radius = 0.4) {
  const lat = THREE.MathUtils.degToRad(latDeg);
  const lon = THREE.MathUtils.degToRad(lonDeg);
  // Earth is offset by earthGroup.position
  const x = Math.cos(lat) * Math.cos(lon) * radius;
  const y = Math.sin(lat) * radius;
  const z = Math.cos(lat) * Math.sin(lon) * radius;
  return new THREE.Vector3(x, y, z).add(earthGroup.position);
}

function showImpactModal(text) {
  const bd = document.getElementById('modalBackdrop');
  const md = document.getElementById('modal');
  const mb = document.getElementById('modalBody');
  if (!bd || !md || !mb) return;
  mb.textContent = text;
  bd.style.display = 'block';
  md.style.display = 'block';
}

window.focusCity = function(lat, lon, name) {
  let marker = cityMarkers.get(name);
  if (!marker) {
    marker = new THREE.Mesh(
      new THREE.SphereGeometry(0.03, 10, 10),
      new THREE.MeshBasicMaterial({ color: 0x66ccff })
    );
    cityMarkers.set(name, marker);
    scene.add(marker);
  }
  marker.position.copy(latLonToVector3(lat, lon));
  // Pan camera target toward Earth and city
  controls.target.copy(earthGroup.position);
  camera.position.lerp(new THREE.Vector3(earthGroup.position.x + 3, 2, earthGroup.position.z + 3), 0.3);
  controls.update();
};

window.showImpactCorridor = function(lat, lon, widthMeters) {
  // Clear old
  while (corridorGroup.children.length) {
    const obj = corridorGroup.children.pop();
    obj.geometry?.dispose?.();
    obj.material?.dispose?.();
  }
  // Draw a simple circle on Earth surface around the city as corridor proxy
  const R = 0.4; // Earth scene radius (scaled mesh radius after scale)
  const angularRadius = Math.min(20, (widthMeters / 1000) / 100) * Math.PI / 180; // crude mapping
  const segments = 128;
  const circlePoints = [];
  const center = latLonToVector3(lat, lon, R);
  // Build circle in tangent plane then project to sphere
  for (let i = 0; i <= segments; i++) {
    const theta = (i / segments) * Math.PI * 2;
    // Local orthonormal basis
    const n = center.clone().sub(earthGroup.position).normalize();
    const t1 = new THREE.Vector3(0,1,0).cross(n).normalize();
    const t2 = n.clone().cross(t1).normalize();
    const pLocal = n.clone().multiplyScalar(R)
      .add(t1.clone().multiplyScalar(Math.cos(theta) * angularRadius * R))
      .add(t2.clone().multiplyScalar(Math.sin(theta) * angularRadius * R));
    circlePoints.push(pLocal);
  }
  const g = new THREE.BufferGeometry().setFromPoints(circlePoints);
  const m = new THREE.LineBasicMaterial({ color: 0xff3366, transparent: true, opacity: 0.6 });
  const loop = new THREE.LineLoop(g, m);
  corridorGroup.add(loop);
};


