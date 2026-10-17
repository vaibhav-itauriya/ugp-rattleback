// 3D view of the rattleback: true semi-ellipsoid on a gridded table, z up.
// Scene units are decimetres (1 unit = 0.1 m) so the camera numbers stay readable.
import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';
import { geom, rotmat } from './model.js';

const U = 10; // metres -> scene units

export class Stage {
  constructor(container) {
    this.el = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.shadowMap.enabled = true;
    this.renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    container.prepend(this.renderer.domElement);
    this.renderer.domElement.setAttribute('aria-label', '3D view of the rattleback; drag to orbit, scroll to zoom');
    this.renderer.domElement.setAttribute('role', 'img');

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(28, 16 / 10, 0.05, 200);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(-4.4, -7.0, 2.3);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.target.set(0, 0, 0.15);
    this.controls.enableDamping = true;
    this.controls.minDistance = 1.2; this.controls.maxDistance = 25;
    this.controls.maxPolarAngle = Math.PI * 0.49;

    this.scene.add(new THREE.HemisphereLight(0xf4f7f6, 0x8a8f8c, 1.1));
    const sun = new THREE.DirectionalLight(0xffffff, 1.7);
    sun.position.set(-3, -2, 8); sun.castShadow = true;
    sun.shadow.mapSize.set(2048, 2048);
    Object.assign(sun.shadow.camera, { left: -4, right: 4, top: 4, bottom: -4, near: 0.5, far: 20 });
    sun.shadow.radius = 4; sun.shadow.bias = -0.0005;
    this.sun = sun; this.scene.add(sun); this.scene.add(sun.target);

    // table: shadow catcher + engineering grid (5 cm minor, 25 cm major)
    const ground = new THREE.Mesh(new THREE.PlaneGeometry(200, 200), new THREE.ShadowMaterial({ opacity: 0.22 }));
    ground.receiveShadow = true; this.scene.add(ground);
    this.grid = new THREE.Group(); this.scene.add(this.grid);
    this.setTheme(document.documentElement.dataset.theme === 'dark' ||
      (!document.documentElement.dataset.theme && matchMedia('(prefers-color-scheme: dark)').matches));

    this.body = new THREE.Group(); this.scene.add(this.body);
    this.marker = new THREE.Mesh(new THREE.RingGeometry(0.035, 0.075, 32), new THREE.MeshBasicMaterial({ color: 0xeb6834, fog: false }));
    this.marker.visible = false; this.scene.add(this.marker);
    this.trail = null; this.data = null; this.exag = 1; this.follow = true; this.showContact = false;
    this.resize = this.resize.bind(this);
    new ResizeObserver(this.resize).observe(container);
    this.resize();
    const loop = () => { this.controls.update(); this.renderer.render(this.scene, this.camera); this.raf = requestAnimationFrame(loop); };
    loop();
  }

  setTheme(dark) {
    this.grid.clear();
    const bg = dark ? 0x10161a : 0xdce3e1;
    this.scene.fog = new THREE.Fog(bg, 7, 16);
    const minor = dark ? 0x1f2a30 : 0xcbd5d2, major = dark ? 0x2e3d44 : 0xaebdb8;
    const mk = (size, div, col) => { const g = new THREE.GridHelper(size, div, col, col); g.rotation.x = Math.PI / 2; g.position.z = 0.0005; return g; };
    this.grid.add(mk(12, 24, minor)); const M = mk(40, 16, major); M.position.z = 0.001; this.grid.add(M);
  }

  resize() {
    const w = this.el.clientWidth, h = this.el.clientHeight;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h; this.camera.updateProjectionMatrix();
  }

  buildBody(p) {
    this.body.clear();
    const wood = new THREE.MeshStandardMaterial({ color: 0x8b5a2b, roughness: 0.55, metalness: 0.0 });
    const cap = new THREE.MeshStandardMaterial({ color: 0xc99566, roughness: 0.7 });
    const ink = new THREE.MeshStandardMaterial({ color: 0x5a3519, roughness: 0.6 });
    const hull = new THREE.SphereGeometry(1, 96, 32, 0, Math.PI * 2, Math.PI / 2, Math.PI / 2);
    hull.rotateX(Math.PI / 2);                    // lower half, opening towards +z
    const mesh = new THREE.Mesh(hull, wood);
    mesh.scale.set(p.a1 * U, p.a2 * U, p.a3 * U); mesh.position.z = p.b * U;
    mesh.castShadow = true;
    const disc = new THREE.Mesh(new THREE.CircleGeometry(1, 96), cap);
    disc.scale.set(p.a1 * U, p.a2 * U, 1); disc.position.z = p.b * U + 1e-4; disc.castShadow = true;
    const stripe = new THREE.Mesh(new THREE.PlaneGeometry(1.84 * p.a1 * U, 0.12 * p.a2 * U), ink);
    stripe.position.z = p.b * U + 2e-3;
    const dot = new THREE.Mesh(new THREE.CircleGeometry(0.22 * p.a2 * U, 32), new THREE.MeshStandardMaterial({ color: 0xf6efe4, roughness: 0.6 }));
    dot.position.set(0.8 * p.a1 * U, 0, p.b * U + 3e-3);           // marks +e1 so the spin sense is visible
    this.body.add(mesh, disc, stripe, dot);
  }

  load(result) {
    this.data = result; const p = result.params;
    this.buildBody(p);
    if (this.trail) { this.scene.remove(this.trail); this.trail.geometry.dispose(); }
    const n = result.T.length, pos = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) { pos[3 * i] = result.Y[9][i] * U; pos[3 * i + 1] = result.Y[10][i] * U; pos[3 * i + 2] = 0.002; }
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    this.trail = new THREE.Line(g, new THREE.LineBasicMaterial({ color: 0x2a78d6 }));
    this.trail.visible = this.showTrail !== false;
    this.scene.add(this.trail);
    this.show(0);
  }

  setTrail(on) { this.showTrail = on; if (this.trail) this.trail.visible = on; }
  setContact(on) { this.showContact = on; this.marker.visible = on && !!this.data; }

  // camera presets, relative to the current target
  setView(name) {
    const t = this.controls.target, off = {
      iso: [-4.4, -7.0, 2.15], front: [0, -8.4, 0.75], side: [-8.4, 0, 0.75], top: [0.001, -0.01, 8.6],
    }[name] || [-4.4, -7.0, 2.15];
    this.camera.position.set(t.x + off[0], t.y + off[1], t.z + off[2]);
    this.controls.update();
  }

  snapshot() { this.renderer.render(this.scene, this.camera); return this.renderer.domElement.toDataURL('image/png'); }

  // place the body at sample i (rocking angles optionally exaggerated about the
  // current heading, so the yaw is untouched and only the tilt is scaled)
  show(i) {
    const R = this.data; if (!R) return;
    const p = R.params, B = { ...p, A2: [p.a1 * p.a1, p.a2 * p.a2, p.a3 * p.a3] };
    let phi = R.Y[3][i], th = R.Y[4][i]; const ps = R.Y[5][i];
    if (this.exag !== 1) { phi *= this.exag; th *= this.exag; }
    const w = [R.Y[0][i], R.Y[1][i], R.Y[2][i]];
    const { h, r } = geom(B, phi, th, ps, w);
    const Q = rotmat(phi, th, ps);
    const m = new THREE.Matrix4().set(
      Q[0][0], Q[0][1], Q[0][2], R.Y[9][i] * U,
      Q[1][0], Q[1][1], Q[1][2], R.Y[10][i] * U,
      Q[2][0], Q[2][1], Q[2][2], h * U,
      0, 0, 0, 1);
    this.body.matrixAutoUpdate = false; this.body.matrix.copy(m); this.body.matrixWorldNeedsUpdate = true;
    if (this.trail) this.trail.geometry.setDrawRange(0, i + 1);
    // computed contact point P = G + Q r, drawn as a ring on the table
    const px = R.Y[9][i] + Q[0][0] * r[0] + Q[0][1] * r[1] + Q[0][2] * r[2];
    const py = R.Y[10][i] + Q[1][0] * r[0] + Q[1][1] * r[1] + Q[1][2] * r[2];
    this.marker.position.set(px * U, py * U, 0.004); this.marker.visible = this.showContact;
    if (this.follow) {
      const tx = R.Y[9][i] * U, ty = R.Y[10][i] * U;
      const dx = tx - this.controls.target.x, dy = ty - this.controls.target.y;
      this.controls.target.x += dx; this.controls.target.y += dy;
      this.camera.position.x += dx; this.camera.position.y += dy;
      this.sun.position.set(tx - 3, ty - 2, 8); this.sun.target.position.set(tx, ty, 0);
      this.grid.position.set(Math.round(tx / 2.5) * 2.5, Math.round(ty / 2.5) * 2.5, 0);   // keep the fine grid under the body
    }
  }
}
