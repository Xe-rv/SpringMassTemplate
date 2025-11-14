using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// Check this out we can require components be on a game object!
[RequireComponent(typeof(MeshFilter))]

public class BParticleSimMesh : MonoBehaviour
{
    public struct BSpring
    {
        public float kd;                        // damping coefficient
        public float ks;                        // spring coefficient
        public float restLength;                // rest length of this spring
        public int attachedParticle;            // index of the attached other particle (use me wisely to avoid doubling springs and sprign calculations)
    }

    public struct BContactSpring
    {
        public float kd;                        // damping coefficient
        public float ks;                        // spring coefficient
        public float restLength;                // rest length of this spring (think about this ... may not even be needed o_0
        public Vector3 attachPoint;             // the attached point on the contact surface
    }

    public struct BParticle
    {
        public Vector3 position;                // position information
        public Vector3 velocity;                // velocity information
        public float mass;                      // mass information
        public BContactSpring contactSpring;    // Special spring for contact forces
        public bool attachedToContact;          // is thi sparticle currently attached to a contact (ground plane contact)
        public List<BSpring> attachedSprings;   // all attached springs, as a list in case we want to modify later fast
        public Vector3 currentForces;           // accumulate forces here on each step        
    }

    public struct BPlane
    {
        public Vector3 position;                // plane position
        public Vector3 normal;                  // plane normal
    }

    public float contactSpringKS = 1000.0f;     // contact spring coefficient with default 1000
    public float contactSpringKD = 20.0f;       // contact spring daming coefficient with default 20

    public float defaultSpringKS = 100.0f;      // default spring coefficient with default 100
    public float defaultSpringKD = 1.0f;        // default spring daming coefficient with default 1

    public bool debugRender = false;            // To render or not to render


    /*** 
     * I've given you all of the above to get you started
     * Here you need to publicly provide the:
     * - the ground plane transform (Transform)
     * - handlePlaneCollisions flag (bool)
     * - particle mass (float)
     * - useGravity flag (bool)
     * - gravity value (Vector3)
     * Here you need to privately provide the:
     * - Mesh (Mesh)
     * - array of particles (BParticle[])
     * - the plane (BPlane)
     ***/

    // Public variables 
    public Transform groundPlaneTransform;      
    public bool handlePlaneCollisions = true;  
    public float particleMass = 1.0f;           
    public bool useGravity = true;            
    public Vector3 gravity = new Vector3(0, -9.8f, 0); 

    // Private variables
    private Mesh mesh;                    
    private BParticle[] particles;           
    private BPlane plane;                    

    /// <summary>
    /// Init everything
    /// HINT: in particular you should probbaly handle the mesh, init all the particles, and the ground plane
    /// HINT 2: I'd for organization sake put the init particles and plane stuff in respective functions
    /// HINT 3: Note that mesh vertices when accessed from the mesh filter are in local coordinates.
    ///         This script will be on the object with the mesh filter, so you can use the functions
    ///         transform.TransformPoint and transform.InverseTransformPoint accordingly 
    ///         (you need to operate on world coordinates, and render in local)
    /// HINT 4: the idea here is to make a mathematical particle object for each vertex in the mesh, then connect
    ///         each particle to every other particle. Be careful not to double your springs! There is a simple
    ///         inner loop approach you can do such that you attached exactly one spring to each particle pair
    ///         on initialization. Then when updating you need to remember a particular trick about the spring forces
    ///         generated between particles. 
    /// </summary>
    void Start()
    {
        // Get the mesh from the mesh filter
        mesh = GetComponent<MeshFilter>().mesh;

        // Initialize particles from mesh vertices
        InitParticles();

        // Initialize the ground plane
        InitPlane();
    }



    /*** BIG HINT: My solution code has as least the following functions
     * InitParticles()
     * InitPlane()
     * UpdateMesh() (remember the hint above regarding global and local coords)
     * ResetParticleForces()
     * ...
     ***/

    void InitParticles()
    {
        Vector3[] vertices = mesh.vertices;
        int vertexCount = vertices.Length;

        particles = new BParticle[vertexCount];

        // Create a particle for each vertex
        for (int i = 0; i < vertexCount; i++)
        {
            particles[i] = new BParticle();
            // Convert from local to world coordinates
            particles[i].position = transform.TransformPoint(vertices[i]);
            particles[i].velocity = Vector3.zero;
            particles[i].mass = particleMass;
            particles[i].attachedToContact = false;
            particles[i].attachedSprings = new List<BSpring>();
            particles[i].currentForces = Vector3.zero;

            // Initialize contact spring with required properties
            particles[i].contactSpring = new BContactSpring
            {
                kd = contactSpringKD,
                ks = contactSpringKS,
                restLength = 0.0f,  
                attachPoint = Vector3.zero
            };
        }

        // Connect particles with springs
        for (int i = 0; i < vertexCount; i++)
        {
            for (int j = i + 1; j < vertexCount; j++)
            {
                // Calculate rest length as initial distance between particles in world space
                float restLen = Vector3.Distance(particles[i].position, particles[j].position);

                // Create spring from i to j
                BSpring spring = new BSpring
                {
                    ks = defaultSpringKS,
                    kd = defaultSpringKD,
                    restLength = restLen,
                    attachedParticle = j
                };

                particles[i].attachedSprings.Add(spring);
            }
        }
    }

    void InitPlane()
    {
        if (groundPlaneTransform != null)
        {
            plane = new BPlane
            {
                position = groundPlaneTransform.position,
                normal = groundPlaneTransform.up.normalized
            };
        }
        else
        {
            // Default ground plane at y=0 facing up
            plane = new BPlane
            {
                position = Vector3.zero,
                normal = Vector3.up
            };
        }
    }

    void ResetParticleForces()
    {
        for (int i = 0; i < particles.Length; i++)
        {
            particles[i].currentForces = Vector3.zero;
        }
    }

    void ApplyGravity()
    {
        if (!useGravity) return;

        for (int i = 0; i < particles.Length; i++)
        {
            particles[i].currentForces += particles[i].mass * gravity;
        }
    }

    void ComputeSpringForces()
    {
        for (int i = 0; i < particles.Length; i++)
        {
            // Iterate through springs attached to this particle
            for (int j = 0; j < particles[i].attachedSprings.Count; j++)
            {
                BSpring spring = particles[i].attachedSprings[j];
                int otherIdx = spring.attachedParticle;

                // Vector from particle i to particle j
                Vector3 x_i = particles[i].position;
                Vector3 x_j = particles[otherIdx].position;
                Vector3 direction = x_j - x_i;
                float currentLength = direction.magnitude;

                // Avoid division by zero
                if (currentLength > 0.0001f) 
                {
                    Vector3 directionNorm = direction / currentLength;

                    // Spring force component: ks * (|xi - xj| - L0) * (xi - xj) / |xi - xj|
                    float displacement = currentLength - spring.restLength;
                    float springMagnitude = spring.ks * displacement;

                    // Damping force component: kd * (vi - vj) · direction
                    Vector3 v_i = particles[i].velocity;
                    Vector3 v_j = particles[otherIdx].velocity;
                    Vector3 relativeVelocity = v_j - v_i;
                    float dampingMagnitude = spring.kd * Vector3.Dot(relativeVelocity, directionNorm);

                    // Total force magnitude
                    float totalMagnitude = springMagnitude + dampingMagnitude;
                    Vector3 force = totalMagnitude * directionNorm;

                    // Apply force to particle i 
                    particles[i].currentForces += force;

                    // Apply reflected force to particle j 
                    particles[otherIdx].currentForces -= force;
                }
            }
        }
    }

    void HandlePlaneCollision()
    {
        if (!handlePlaneCollisions) return;

        // Update plane position and normal from transform
        if (groundPlaneTransform != null)
        {
            plane.position = groundPlaneTransform.position;
            plane.normal = groundPlaneTransform.up.normalized;
        }

        for (int i = 0; i < particles.Length; i++)
        {
            // Calculate signed distance from plane
            Vector3 toParticle = particles[i].position - plane.position;
            float distance = Vector3.Dot(toParticle, plane.normal);

            // Check if penetrating (distance < 0 means particle is on negative side of plane)
            if (distance < 0)
            {
                // If not already attached, initialize the contact spring
                if (!particles[i].attachedToContact)
                {
                    particles[i].attachedToContact = true;

                    // Find nearest point on plane at moment of contact
                    // This is the attach point and should remain fixed during penetration
                    particles[i].contactSpring.attachPoint = particles[i].position - distance * plane.normal;

                    // Set contact spring properties (ks=1000, kd=20, restLength=0)
                    particles[i].contactSpring.ks = contactSpringKS;
                    particles[i].contactSpring.kd = contactSpringKD;
                    particles[i].contactSpring.restLength = 0.0f;
                }

                // Compute penalty force: -ks * (xp - xg) · n * n - kd * vp
                Vector3 x_p = particles[i].position;
                Vector3 x_g = particles[i].contactSpring.attachPoint;
                Vector3 v_p = particles[i].velocity;

                // Spring component: -ks * (xp - xg) · n * n
                Vector3 penetrationVector = x_p - x_g;
                float penetrationDepth = Vector3.Dot(penetrationVector, plane.normal);
                Vector3 springForce = -contactSpringKS * penetrationDepth * plane.normal;

                // Damping component: -kd * vp (full velocity, not just normal component based on image)
                Vector3 dampingForce = -contactSpringKD * v_p;

                particles[i].currentForces += springForce + dampingForce;
            }
            else
            {
                // No longer penetrating, detach contact spring
                if (particles[i].attachedToContact)
                {
                    particles[i].attachedToContact = false;
                }
            }
        }
    }

    void IntegrateParticles(float dt)
    {
        for (int i = 0; i < particles.Length; i++)
        {
            // Calculate acceleration from forces: a = F / m
            Vector3 acceleration = particles[i].currentForces / particles[i].mass;

            // Symplectic Euler: update velocity first using current forces
            particles[i].velocity += acceleration * dt;

            // Then update position using the NEW velocity
            particles[i].position += particles[i].velocity * dt;
        }
    }

    void UpdateMesh()
    {
        Vector3[] vertices = mesh.vertices;

        for (int i = 0; i < particles.Length; i++)
        {
            // Convert from world to local coordinates
            vertices[i] = transform.InverseTransformPoint(particles[i].position);
        }

        mesh.vertices = vertices;

        // Update bounds and normals as required
        mesh.RecalculateBounds();
        mesh.RecalculateNormals();
    }

    void FixedUpdate()
    {
        float dt = Time.fixedDeltaTime;

        // Reset forces
        ResetParticleForces();

        // Apply all forces
        ApplyGravity();
        ComputeSpringForces();
        HandlePlaneCollision();

        // Integrate
        IntegrateParticles(dt);

        // Update mesh
        UpdateMesh();
    }

    /// <summary>
    /// Draw a frame with some helper debug render code
    /// </summary>
    public void Update()
    {
        if (debugRender)
        {
            int particleCount = particles.Length;
            for (int i = 0; i < particleCount; i++)
            {
                Debug.DrawLine(particles[i].position, particles[i].position + particles[i].currentForces, Color.blue);

                int springCount = particles[i].attachedSprings.Count;
                for (int j = 0; j < springCount; j++)
                {
                    Debug.DrawLine(particles[i].position, particles[particles[i].attachedSprings[j].attachedParticle].position, Color.red);
                }
            }
        }
       
    }
}