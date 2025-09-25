import * as THREE from "three";

export default function getStarfield({ numStars = 1000, radius = 200 } = {}) {
  const geometry = new THREE.BufferGeometry();
  const positions = new Float32Array(numStars * 3);
  for (let i = 0; i < numStars; i++) {
    // random point on a sphere shell
    const theta = Math.acos(THREE.MathUtils.randFloatSpread(2));
    const phi = Math.random() * Math.PI * 2;
    const r = radius;
    const x = r * Math.sin(theta) * Math.cos(phi);
    const y = r * Math.cos(theta);
    const z = r * Math.sin(theta) * Math.sin(phi);
    positions.set([x, y, z], i * 3);
  }
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  const material = new THREE.PointsMaterial({ color: 0xffffff, size: 0.7, sizeAttenuation: true, transparent: true, opacity: 0.9 });
  return new THREE.Points(geometry, material);
}


